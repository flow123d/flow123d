import os
import xml.etree.ElementTree as ET
import argparse

def create_junit_xml(status, xml_output, log, class_name):
    if status == 0:
        return

    os.makedirs(os.path.dirname(xml_output), exist_ok=True)

    testsuites = ET.Element("testsuites")
    testsuite = ET.SubElement(
        testsuites, "testsuite", 
        name=class_name, 
        tests="1", 
        failures="1" if status != 0 else "0",
        errors="0", 
        time="0.0"
    )
    testcase = ET.SubElement(
        testsuite, "testcase", 
        classname=class_name, 
        name=f"compilation-failed-for-{class_name}", 
        time="0.0"
    )

    failure = ET.SubElement(
        testcase, "failure", 
        message=f"Compilation failed for {class_name}", 
        type=f"compilation-failed-for-{class_name}",
        text=log
    )
    failure.text = log

    tree = ET.ElementTree(testsuites)
    with open(xml_output, "wb") as xml_file:
        tree.write(xml_file, encoding="utf-8", xml_declaration=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate JUnit XML report for compilation")
    parser.add_argument("--status", type=int, required=True, help="Compilation status")
    parser.add_argument("--output", type=str, required=True, help="Path to the output XML file")
    parser.add_argument("--log", type=str, required=True, help="Compilation output and error log")
    parser.add_argument("--class-name", type=str, default="Compilation", help="Class name for the test case")
    
    args = parser.parse_args()
    create_junit_xml(args.status, args.output, args.log, args.class_name)
