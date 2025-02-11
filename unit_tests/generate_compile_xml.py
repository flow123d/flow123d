import sys
import os
import xml.etree.ElementTree as ET

def create_junit_xml(status, log_file, xml_output):
    with open(log_file, "r", encoding="utf-8") as f:
        compile_output = f.read()

    os.makedirs(os.path.dirname(xml_output), exist_ok=True)

    testsuites = ET.Element("testsuites")
    testsuite = ET.SubElement(testsuites, "testsuite", name="CompilationTests", tests="1", failures="1" if status != 0 else "0", errors="0", time="0.0")
    testcase = ET.SubElement(testsuite, "testcase", classname="Compilation", name="CompileTest", time="0.0")

    if status == 127:
        failure = ET.SubElement(testcase, "failure", message="Test binary missing", type="RuntimeError")
        failure.text = compile_output
    elif status != 0:
        failure = ET.SubElement(testcase, "failure", message="Compilation failed", type="CompileError")
        failure.text = compile_output
    else:
        system_out = ET.SubElement(testcase, "system-out")
        system_out.text = compile_output

    tree = ET.ElementTree(testsuites)
    with open(xml_output, "wb") as xml_file:
        tree.write(xml_file, encoding="utf-8", xml_declaration=True)

if __name__ == "__main__":
    status = int(sys.argv[1])
    log_file = sys.argv[2]
    xml_output = sys.argv[3]

    create_junit_xml(status, log_file, xml_output)
    sys.exit(status)