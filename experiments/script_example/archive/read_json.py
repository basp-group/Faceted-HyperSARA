import json

if __name__ == "__main__":

    with open("experiments/script_example/setup_test.json") as f:
        data = json.load(f)

    print(data)
