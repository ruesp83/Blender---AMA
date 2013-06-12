import os
import sys
import string


def _main(argv):
    if len(sys.argv) < 2:
        sys.stdout.write("Usage: create-co <data_file>\n")
        sys.exit(1)

    filename = sys.argv[1]

    try:
        fpin = open(filename, "r")
    except:
        sys.stdout.write("Unable to open input %s\n" % sys.argv[1])
        sys.exit(1)

    try:
        fpout = open("copy_patch.bat", "w")
    except:
        sys.stdout.write("Unable to open output %s\n" % cname)
        sys.exit(1)

    for line in iter(fpin):
        if (string.find(line, "Index: ", 0, 7) != -1):
            path = string.replace(line[7:], "/", "\\")
            path = path[0:len(path)-1]
            new_line = "copy " + path + " J:\\Blender---AMA\\"
            path = path[0:string.rfind(path, "\\")]
            new_line = new_line + path + " /Y"
            fpout.write(new_line + "\n")
    fpout.write("copy patch-ama J:\Blender---AMA\patch-ama /Y\n")
    sys.stdout.write("Files Created: copy_patch.bat\n")
    fpout.close()
    fpin.close()


if __name__ == "__main__":
    _main(sys.argv)