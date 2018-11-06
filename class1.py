# Interpreter: Python 3.7
# Author: Risheek Rajolu
# Date: August 28, 2018

import numpy

myList = ["apple", "banana"]
myList.append("orange");
myList.extend(["pineapple", "coconut"])
myList.insert(0, "artichoke")
print(myList)

myDict = {"bob":31, "joan":29}
myDict["children"] = {"mary": 3, "gary": 4}
print(list(myDict.keys()))

array = numpy.array([2,3,4])
print("mean: " + str(array.mean()))
print("std: " + str(array.std()))

print("myDict: {")
for x in myDict.keys():
    if isinstance(myDict[x], dict):
        print("\t" + str(x) + " - {")
        for y in myDict[x]:
            print("\t\t" + str(y) + " - " + str(myDict[x][y]))
        print("\t}")
    else:
        print("\t" + str(x) + " - " + str(myDict[x]))
print("}")
