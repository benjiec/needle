import requests
import re


def fetch_module_data(module_id):
    url = f"https://rest.kegg.jp/get/{module_id}"
    module_data = requests.get(url).text
    return module_data


def parse_module_definition(module_data):
    keyword = "DEFINITION"
    lines = module_data.split("\n")
    def_lines = [line for line in lines if line.startswith(keyword)]
    if len(def_lines) != 1:
        return None
    def_line = def_lines[0][len(keyword):]
    return [x for x in re.split(r"[\(\)\W\s\,]+", def_line) if x]


module_ko_numbers = {}

with open('data/modules.tsv') as f:
    for line in f.readlines():
        if line.startswith("Module ID"):
            continue
        module_id, module_name = line.split("\t")
        module_data = fetch_module_data(module_id)
        ko_numbers = parse_module_definition(module_data)
        module_ko_numbers[module_id] = ko_numbers
        print(module_id)

with open('data/module_ko.tsv', "w+") as f:
    f.write("Module ID\tModule Definition\n")
    for module_id, ko_numbers in module_ko_numbers.items():
        ko_number_str = ','.join(ko_numbers)
        f.write(f"{module_id}\t{ko_number_str}\n")
