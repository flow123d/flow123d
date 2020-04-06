#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import sys
import os
import re
from typing import List, Dict


is_debug = str(os.environ.get("DEBUG", "0")) in ("1", "true")

class Col:
    def __init__(self, name: str, cols: int, func: any):
        self.name = name
        self.cols = cols
        self.func = func
        self.children = []
        self.fmt = "{result:%ds}" % self.cols
        self.fmt_header = "{name:^%ds}" % (self.cols)

    def __call__(self, context: Dict):
        result = str(self.func(context))
        return self.fmt.format(result=result)
    
    def header(self):
        return self.fmt_header.format(name=self.name)

def extract_tag(context: Dict):
    return "{indent}{tag}".format(**context)

def extract_dur(context: Dict):
    try:
        time = float(context["cumul-time-sum"] / context["call-count-sum"])
    except:
        time = 0.0

    return "{time:7,d}".format(time=int(time*1000)).replace(',', ' ')

def extract_calls(context: Dict):
    return "{call-count-sum:5d}".format(**context)

def extract_rel_per(context: Dict):
    return "{percent:3.0f}".format(**context)

def extract_abs_per(context: Dict):
    return "{abs_per:3.0f}".format(**context)

def extract_ch_tot(context: Dict):
    return "{children_total:6d}".format(**context)

def extract_src(context: Dict):
    file_path = str(context.get("file-path", ""))
    file_path = file_path[5:] if file_path.startswith("/src/") else file_path
    file_line = context.get("file-line")
    file_path = './' if not file_path else file_path
    if file_line:
        return "{file_path}:{file_line}".format(file_path=file_path, file_line=context.get("file-line"))
    return "{file_path}".format(file_path=file_path)

def extract_function(context: Dict):
    return "{function}".format(**context)

def extract_min_max(context: Dict):
    try:
        ratio = float(context.get("cumul-time-min") / context.get("cumul-time-max"))
    except:
        ratio = 0.0
    return "{ratio:6.3f}".format(ratio=ratio)

def extract_rem(context: Dict):
    if not len(context["children"]):
        return"{remainder:6s}".format(remainder="")

    remainder = context["time"] - context["children_total"]
    return "{remainder:6d}".format(remainder=remainder)

def extract_rem_prc(context: Dict):
    if not len(context["children"]):
        return"{remainder:3s}".format(remainder="")

    remainder = (context["children_total"] / context["safe_time"]) * 100
    return "{remainder:3.0f}".format(remainder=remainder)

def enrich_nodes(node: Dict, parent: Dict, root: Dict, indent=0):
    node["children"] = node.get("children", [])
    node["len_children"] = str(len(node["children"] ) or " ")
    node["time"] = int(node.get("cumul-time-sum") * 1000)
    node["safe_time"] = max(node["time"], 1)
    node["indent"] = indent * ' '
    node["level"] = indent
    node["parent_time"] = parent["time"]
    node["abs_per"] = node["time"] * 100 / root["safe_time"]

    for idx, subnode in enumerate(node["children"]):
        if indent < 15:
            enrich_nodes(subnode, node, root, indent+1)

def longest_common_prefix(strs):
    longest_pre = ""
    if not strs: return longest_pre
    shortest_str = min(strs, key=len)
    for i in range(len(shortest_str)):
        if all([x.startswith(shortest_str[:i+1]) for x in strs]):
            longest_pre = shortest_str[:i+1]
        else:
            break
    return longest_pre

def squashed_node(rest: List[Dict]):
    if len(rest) == 1:
        return rest[0]

    def sum_prop(p: str, f=sum):
        return f([x[p] for x in rest])

    base = rest[0].copy()
    fps = list(sum_prop("file-path", set))
    fnc = list(sum_prop("function", set))
    base.update(**{
        "tag": "others (%d more)" % len(rest),
        "time": sum_prop("time"),
        "call-count-sum": sum_prop("call-count-sum"),
        "percent": sum_prop("percent"),
        "abs_per": sum_prop("abs_per"),
        "function": longest_common_prefix(fnc), # try to get single result if possible
        "file-path": longest_common_prefix(fps), # try to get single result if possible
        "file-line": "",
    })
    return base

def traverse_children(node: Dict, filter: any, sort_function: any):
    children = sort_function(node["children"])
    rest = []
    for subnode in children:
        if filter(subnode):
            yield subnode
        else:
            rest.append(subnode)
    
    if rest:
        yield squashed_node(rest)

def visit_node(node: Dict, parent: Dict, root: Dict, filter: any, sort_function, columns: List[Col], indent=1):
    node["children_total"] = sum([n["time"] for n in node["children"]])

    lines = [' | '.join([col(node) for col in columns])]
    for idx, subnode in enumerate(traverse_children(node, filter, sort_function)):
        if indent < 15:
            lines.extend(visit_node(subnode, node, root, filter, sort_function, columns, indent+1))
    return lines

def get_env(feature: str, default=None):
    key = "FLOW123D_PROFILER_{feature}".format(feature=feature.upper())
    return os.environ.get(key, default)

def is_on(feature: str, default=True):
    if is_debug:
        return True
    key = "FLOW123D_PROFILER_{feature}".format(feature=feature.upper())
    value = str(os.environ.get(key, default)).lower()
    return value in ("1", "true")


class SimpleTableFormatter2 (object):
    """
    Class which takes json object from flow123d benchmark profiler report
     and returns simple table-like text string
    """

    def __init__(self):
        columns = [
            Col('Tag', 55, extract_tag) if is_on("TAG") else None,
            Col('∑ Time', 7, extract_dur) if is_on("TIME") else None,
            Col('%Cov', 4, extract_rem_prc) if is_on("COV") else None,
            Col('# Cls', 5, extract_calls) if is_on("CLS") else None,
            Col('min/max', 7, extract_min_max) if is_on("REL") else None,
            Col('%Rel', 4, extract_rel_per) if is_on("REL") else None,
            Col('%Abs', 4, extract_abs_per) if is_on("ABS") else None,
            Col('Func', 28, extract_function) if is_on("FUNC", False) else None, # env var FLOW123D_PROFILER_FUNC must be 1
            Col('Src', 40, extract_src) if is_on("SRC") else None,
        ]
        self.columns = [c for c in columns if c]

    def set_styles(self, styles):
        """Overrides default styles"""
        pass

    def _get_header(self, json: Dict, title: str, value=None):
        value = value or json.get(title.lower().replace(' ', '-'))
        return "{title:20s} {value}".format(title=title, value=value)

    def get_summary(self, json: Dict):
        yield self._get_header(json, "Program name", "{program-name} ({program-version})".format(**json))
        yield self._get_header(json, "Program branch", "{program-branch} ({program-revision})".format(**json))
        yield self._get_header(json, "Program build")

        if 'source-dir' in json:
            yield self._get_header(json, "Source dir")

        yield self._get_header(json, "Timer resolution", "{timer-resolution:1.6f}".format(**json))

        desc = re.sub("\s+", " ", json["task-description"], re.M)
        yield self._get_header(json, "Task description", desc)
        yield self._get_header(json, "Task size")

        yield self._get_header(json, "Run process count")
        yield self._get_header(json, "Run started", json["run-started-at"])
        yield self._get_header(json, "Run ended", json["run-finished-at"])

    def format(self, json):
        """"Formats given json object"""
        lines = []

        # determine squashin condition via FLOW123D_PROFILER_SQUASH_OFF
        squash_off = is_on('SQUASH_OFF', False)
        if squash_off:
            is_significant = lambda x: True
        else:
            is_significant = lambda x: x["abs_per"] > 1.0 or x["time"] > 100

        # determine sort field via FLOW123D_PROFILER_SORT
        sort_field = get_env("SORT", "cumul-time-sum")
        sort_function = lambda x: sorted(x, key=lambda y: y.get(sort_field, 0), reverse=True)

        # determine summary via FLOW123D_PROFILER_SUMMARY
        # and possible append it to the output
        header = [col.header() for col in self.columns]
        if is_on("SUMMARY"):
            lines.extend(self.get_summary(json))
            lines.append("")
            lines.append('-x-'.join([len(h) * '-' for  h in header]))

        # add table start
        lines.extend([
            ' | '.join(header),
            '-x-'.join([len(h) * '-' for  h in header])
        ])

        # get the first node
        first = json['children'][0]
        # add some other metrics
        enrich_nodes(first, first, first)

        # append the rest of the table
        lines.extend(
            visit_node(
                first, first, first,
                filter=is_significant,
                sort_function=sort_function,
                columns=self.columns,
            )
        )

        # join by newline and viola
        return '\n'.join(lines)
