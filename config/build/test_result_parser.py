import argparse
import os
import glob
import xml.etree.ElementTree as ET

# python .\xml_process.py test-results fem 12284539523 34280995371

# Example of URL for test row:
# https://github.com/MocStepan/flow123d/actions/runs/12282944892/job/34275773206#step:6:7
# step means which step in the job
# after step there is row number
# runs/{run_id}/job/{job_id}#step:{step}:{row}
# I will need to extract row from logs but firstly i need to regex the logs
# based on how much cores is used for tests. Core number is in the xml test name
# Xml test name is formatted as: test-name-{core_number}.xml
# Min cores can be 1 and max can be 3
# After i get the core number i take the regex group based on the core number 
# and extract the row number
# then i will use row number to create url and after that i will add url to xml file
# If xml contain failed test.
# It will create new directory for only failed tests.

TEST_STEP_NUMBER = 5
GITHUB_URL = "https://github.com/MocStepan/flow123d/actions/runs/"

def get_url(run_id, job_id, row):
    return f"{GITHUB_URL}{run_id}/job/{job_id}#step:{TEST_STEP_NUMBER}:{row}"

def get_row_number(log_file, core_number, test_name):
    line_number = 1
    for index, line in enumerate(log_file.splitlines(), start=1):
        if line.startswith(f"[ RUN      ] {test_name}"):
            if line_number == core_number:
                return index
            else:
                line_number += 1

    print(f"Test '{test_name}' not found in log.")
    return None


def process(xml_file, log_file, run_id, job_id):
    core_number = int(xml_file.split("-")[-1].split(".")[0])
    tree = ET.parse(xml_file)
    root = tree.getroot()

    for test_suite in root.findall("testsuite"):
        for test_case in test_suite.findall("testcase"):
            if test_case.find("failure") is not None:
                test_name = f"{test_case.get('classname')}.{test_case.get('name')}"
                row_number = get_row_number(log_file, core_number, test_name)
                url = get_url(run_id, job_id, row_number)
                print(f"Test '{test_name}' failed. URL: {url}")
    return



def add_url_to_xml(xml_file, url):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    return



def process_xml_files(build_dir, test_dir, run_id, job_id):
    print(f"Processing XML files for Run ID: {run_id}, Job ID: {job_id}")
    print(f"Build directory: {build_dir}, Test directory: {test_dir}")
    
    xml_path = os.path.join(build_dir, "unit_tests", test_dir, "*.xml")
    xml_files = glob.glob(xml_path)

    log_file_path = os.path.join(build_dir, "unit_tests", test_dir, f"{test_dir}_log.txt")
    with open(log_file_path, 'r', encoding='utf-8') as log_file:
        log_content = log_file.read()

    if not xml_files:
        print(f"No XML files found in directory: {xml_path}")
        return

    for xml_file in xml_files:
        print(f"Processing file: {xml_file}")

        try:
            process(xml_file, log_content, run_id, job_id)
        except ET.ParseError as e:
            print(f"Failed to parse {xml_file}: {e}")

   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process XML files from test directories.")
    parser.add_argument("build_dir", help="Path to the build directory")
    parser.add_argument("test_dir", help="Test directory to process XML files from")
    parser.add_argument("run_id", help="GitHub Actions Run ID")
    parser.add_argument("job_id", help="GitHub Actions Job ID")

    args = parser.parse_args()

    process_xml_files(args.build_dir, args.test_dir, args.run_id, args.job_id)