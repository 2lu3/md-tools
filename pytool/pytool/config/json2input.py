import json

def json2input(json_text: str) -> str:
    result: str = ""
    for section, values in json.loads(json_text).items():
        result += " " * 8
        result += f'f.write("[{section}]\\n")\n'
        for key, value in values.items():
            result += " " * 8
            if "{" in value:
                result += f'f.write(f"{key:<20} = {value}\\n")\n'
            else:
                result += f'f.write("{key:<20} = {value}\\n")\n'

        result += "\n"
    return result



