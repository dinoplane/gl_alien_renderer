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

        DIST_INCR = 2
        ANGLE_INCR = 5
        SCALE = 1.0

        for i in range(width):
            for j in range(height):
                for k in range(depth):
                    entData = {
                        "classname": "peeper",
                        "origin": "{} {} {}".format(DIST_INCR*i, DIST_INCR*j, DIST_INCR*k),
                        "angles": "{} {} {}".format(ANGLE_INCR*i, ANGLE_INCR*j, ANGLE_INCR*k),
                        
                        "scale": "{} {} {}".format(SCALE, SCALE, SCALE),
                        
                        "mesh" : "./resources/assets/models/arctic_peeper/scene.gltf",
                        "material": "base_inst.shader",
                        "is_instanced" : "1"
                    }
                    dump_entity_data(f, entData)

        dump_entity_data(f, {
            "classname" : "camera",
            "origin" : "-3 3 -3",
            "angles" : "-60 -315 0",
            "fovy" : "30",
            "near" : "0.1",
            "far" : "10.0",
        })

        dump_entity_data(f, {
            "classname" : "camera",
            "origin" : "-10 20 -10",
            "angles" : "-60 -315 0",
            "fovy" : "30",
            "near" : "0.1",
            "far" : "100.0",
        })


if __name__ == '__main__':
    main()