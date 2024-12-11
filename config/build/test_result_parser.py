import argparse
import os
import glob
import xml.etree.ElementTree as ET

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
            # Example: Extract and print all tags and text
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