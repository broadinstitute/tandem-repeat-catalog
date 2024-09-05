"""Basic validation for any JSON file. Makes sure that the file is in valid JSON format and that it contains a list
of dictionaries that have the required keys."""

import argparse
import ijson
import gzip
import os

def failed_validation(json_path, keys=None):
	keys = set(keys) if keys is not None else set()
	error_counter =	total = 0
	fopen = gzip.open if json_path.endswith("gz") else open
	with fopen(json_path, "rb") as f:
		parser = ijson.items(f, "item", use_float=True)
		for i, record in enumerate(parser):
			total += 1
			if not isinstance(record, dict):
				print(f"ERROR: record #{i + 1} is not a dictionary: {record}")
				error_counter += 1
				if error_counter >= 50:
					return error_counter, total

			if keys:
				missing_keys = keys - set(record.keys())
				if missing_keys:
					print(f"ERROR: {', '.join(keys)} key(s) are missing in record #{i + 1}: {record}")
					error_counter += 1
					if error_counter >= 50:
						return error_counter, total

	return error_counter, total


def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
	parser.add_argument("-k", "--key", action="append", help="Key to check for in each record")
	parser.add_argument("json_path", help="Path of the JSON file to validate")
	args = parser.parse_args()

	if not os.path.isfile(args.json_path):
		parser.error(f"File not found: {args.json_path}")

	error_counter, total = failed_validation(args.json_path, keys=args.key)
	if error_counter == 0:
		print(f"All {total:,d} rows PASSED validation")
	else:
		print(f"ERROR: {error_counter:,d} out of {total:,d} ({error_counter/total:.1%}) rows FAILED validation")
		exit(1)


if __name__ == "__main__":
	main()