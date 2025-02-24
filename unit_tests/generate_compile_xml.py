import os
import xml.etree.ElementTree as ET
import argparse

def create_junit_xml(status, xml_output, url, repository, sha, run_id):
    if status == 0:
        return

    os.makedirs(os.path.dirname(xml_output), exist_ok=True)

    testsuites = ET.Element("testsuites")
    testsuite = ET.SubElement(testsuites, "testsuite", name="CompilationTests", tests="1", failures="1" if status != 0 else "0", errors="0", time="0.0")
    testcase = ET.SubElement(testsuite, "testcase", classname="Compilation", name="CompileTest", time="0.0")
    failure = ET.SubElement(testcase, "failure", message="Compilation failed", type="CompileError")
    failure.text = f"{url}/{repository}/commit/{sha}/checks/{run_id}"

    tree = ET.ElementTree(testsuites)
    with open(xml_output, "wb") as xml_file:
        tree.write(xml_file, encoding="utf-8", xml_declaration=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate JUnit XML report")
    parser.add_argument("--status", type=int, required=True, help="Test status")
    parser.add_argument("--output", type=str, required=True, help="Path to the output XML file")
    parser.add_argument("--url", type=str, required=True, help="URL of the check")
    parser.add_argument("--repository", type=str, required=True, help="Repository name")
    parser.add_argument("--sha", type=str, required=True, help="Commit SHA")
    parser.add_argument("--run_id", type=str, required=True, help="Run ID")
    
    args = parser.parse_args()
    create_junit_xml(args.status, args.output, args.url, args.repository, args.sha, args.run_id)