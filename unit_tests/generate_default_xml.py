#!/usr/bin/env python3
import sys
import os
import xml.etree.ElementTree as ET

def create_junit_xml(compile_output, compile_errors, xml_path):
    testsuites = ET.Element("testsuites")
    
    testsuite_attrs = {
        "name": "Compilation",
        "tests": "1",
        "failures": "1" if compile_errors.strip() else "0",
        "errors": "0",
        "time": "0"
    }
    testsuite = ET.SubElement(testsuites, "testsuite", testsuite_attrs)
    
    testcase = ET.SubElement(testsuite, "testcase", {
        "classname": "Compilation",
        "name": "CompileTest"
    })
    
    if compile_errors.strip():
        failure = ET.SubElement(testcase, "failure", {
            "message": "Compilation errors occurred",
            "type": "CompilationError"
        })
        failure.text = compile_errors.strip()
    else:
        system_out = ET.SubElement(testcase, "system-out")
        system_out.text = compile_output.strip()

    tree = ET.ElementTree(testsuites)
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)
    print(f"XML report byl úspěšně vygenerován: {xml_path}")

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: generate_default_xml.py <compile_output.txt> <compile_errors.txt> <output_report.xml>")
    
    compile_output_file = sys.argv[1]
    compile_errors_file = sys.argv[2]
    xml_report_file = sys.argv[3]

    compile_output = ""
    compile_errors = ""
    if os.path.exists(compile_output_file):
        with open(compile_output_file, "r") as f:
            compile_output = f.read()
    if os.path.exists(compile_errors_file):
        with open(compile_errors_file, "r") as f:
            compile_errors = f.read()

    create_junit_xml(compile_output, compile_errors, xml_report_file)

if __name__ == '__main__':
    main()
