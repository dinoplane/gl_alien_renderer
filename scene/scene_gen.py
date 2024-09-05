import sys

def dump_entity_data(file_handle, entity_data):
    file_handle.write("{\n")
    for key, value in entity_data.items():
        file_handle.write("\"{}\" \"{}\"\n".format(key, value))
    file_handle.write("}\n")

def main():
    print("Hello World")
    print(sys.argv)

    with open(sys.argv[1], 'w') as f:
        width = int(sys.argv[2])
        height = int(sys.argv[3])
        depth = int(sys.argv[4])

        for i in range(width):
            for j in range(height):
                for k in range(depth):
                    entData = {
                        "classname": "cube",
                        "origin": "{} {} {}".format(2*i, 2*j, 2*k)
                    }
                    dump_entity_data(f, entData)

        dump_entity_data(f, {
            "classname" : "camera",
            "origin" : "-3 3 -3",
            "angles" : "-315 -60 0",
            "fovY" : "30",
            "near" : "0.1",
            "far" : "10.0",
        })

        dump_entity_data(f, {
            "classname" : "camera",
            "origin" : "-10 20 -10",
            "angles" : "-315 -60 0",
            "fovY" : "30",
            "near" : "0.1",
            "far" : "100.0",
        })


if __name__ == '__main__':
    main()