import argparse
import os
import glob
import xml.etree.ElementTree as ET

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

def get_core_number(test_name):
    return

def get_row_number(logs, core_number):
    return

def get_url(run_id, job_id, step, row):
    return

def add_url_to_xml(xml_file, url):
    return


def process_xml_files(build_dir, test_dir, job_id, run_id):
    print(f"Processing XML files for Job ID: {job_id}, Run ID: {run_id}")
    print(f"Build directory: {build_dir}, Test directory: {test_dir}")
    xml_path = os.path.join(build_dir, "unit_tests", test_dir, "*.xml")
    xml_files = glob.glob(xml_path)

    if not xml_files:
        print(f"No XML files found in directory: {xml_path}")
        return

    print(f"Processing XML files for Job ID: {job_id}, Run ID: {run_id}")

    for xml_file in xml_files:
        print(f"Processing file: {xml_file}")
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            for element in root.iter():
                print(f"Tag: {element.tag}, Text: {element.text}")
        except ET.ParseError as e:
            print(f"Failed to parse {xml_file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process XML files from test directories.")
    parser.add_argument("build_dir", help="Path to the build directory")
    parser.add_argument("test_dir", help="Test directory to process XML files from")
    parser.add_argument("job_id", help="GitHub Actions Job ID")
    parser.add_argument("run_id", help="GitHub Actions Run ID")

    args = parser.parse_args()

    process_xml_files(args.build_dir, args.test_dir, args.job_id, args.run_id)