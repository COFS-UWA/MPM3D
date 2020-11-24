import os

#file1_name = "t3d_mt_model.txt"
#file2_name = "t3d_s_model.txt"
file1_name = "t3d_stp_mt.txt"
file2_name = "t3d_mt_model.txt"

file1 = open("../Build/TestsParallel/" + file1_name, "r")
file2 = open("../Build/TestsParallel/" + file2_name, "r")

line_num = 0
is_the_same = True
line1 = file1.readline()
line2 = file2.readline()
while True:
    line_num += 1
    if line1 != line2:
        is_the_same = False
        break
    if not line1 or not line2:
        break
    line1 = file1.readline()
    line2 = file2.readline()

print("%s and %s are " % (file1_name, file2_name), end="")
if not is_the_same:
    print("different at line %d." % line_num)
else:
    print("the same.")

file1.close();
file2.close();
