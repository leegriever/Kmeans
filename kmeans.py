import math
import sys


def main():
    try:
        K = int(sys.argv[1])
    except:
        printInvalidInput()

    numOfArgs = len(sys.argv)
    if numOfArgs == 4:
        max_iter = 200
        input = sys.argv[2]
        output = sys.argv[3]
    elif numOfArgs == 5:
        try:
            max_iter = int(sys.argv[2])
        except:
            printInvalidInput()
        input = sys.argv[3]
        output = sys.argv[4]
    else:
        printInvalidInput()
    
    try:
        vectors = readInput(input)

    except:
        printInvalidInput()
    
    validity(K, len(vectors), max_iter)
    
    try:
        centroids = intializeCentroids(vectors,K)

        for i in range(max_iter):
            clusters = [[] for i in range(K)]
            num_of_unchanged_centorids = 0
            for j in range(len(vectors)):
                assignVectorToCluster(vectors,centroids,clusters,j)
            for l in range(K):
                old_centroid = centroids[l].copy()
                updateCentroid(vectors,centroids,clusters,l)
                new_centroid = centroids[l]
                if math.sqrt(Distance(old_centroid,new_centroid))<0.001:
                    num_of_unchanged_centorids += 1
            if num_of_unchanged_centorids==K:
                break

        writeToOutput(output,centroids)  
    except:
        print("An Error Has Occurred")
        sys.exit(1)
          
            

def printInvalidInput():
    print("Invalid Input!")
    sys.exit(1)


def validity(k, n, max_iter):
    if k <= 1 or k >= n:
        printInvalidInput()
    if (k % 1) != 0:
        printInvalidInput()
    if n <= 0:
        printInvalidInput()
    if (max_iter % 1) != 0:
        printInvalidInput()
    if max_iter <= 0:
        printInvalidInput()
    

# func to read data from file, return a list of lstsc(vectors = [],vectorCooirdinates = [])
def readInput(fileName):
    file = open(fileName,'r')
    lines = file.read().splitlines()
    file.close
    vectors = []
    for i in range(len(lines)):
        line = lines[i].split(',')
        v = []
        for j in range(len(line)):
            coordinate = float(line[j])
            v.append(coordinate)
        vectors.append(v)

    return vectors  


# creates lst in length k, initialize with first k vectors
def intializeCentroids(vectors,k):
    centroids = []
    for i in range(k):
        copy_vector = (vectors[i]).copy()
        centroids.append(copy_vector)
    return centroids

# calculates euclideanDistance between 2 vectors, returns dis
def Distance(x,y):
    sumOfSquaredDiffrence = 0
    for i in range(len(x)):
        sumOfSquaredDiffrence += math.pow(float(x[i])-float(y[i]),2)
    return sumOfSquaredDiffrence


# updates and returns centroid
def updateCentroid(vectors,centroids,clusters,i):
    cluster = clusters[i]
    cluster_size = len(cluster)
    for j in range(len(vectors[0])):
        sum = 0.0
        for vectorIndex in cluster:
            sum += vectors[vectorIndex][j]
        centroids[i][j] = sum/cluster_size
    return


# assigns vector to closest cluster
def assignVectorToCluster(vectors,centroids,clusters,vectorIndex):
    min = sys.float_info.max
    index = 0
    for i in range (len(centroids)):
        distance = Distance(vectors[vectorIndex],centroids[i])
        if distance<min:
            min = distance
            index = i
    clusters[index].append(vectorIndex)
    return

def writeToOutput(fileName,centroids):
    file = open(fileName,'w')
    num_of_coordinates_in_vector = len(centroids[0])
    for i in range(len(centroids)):
        for j in range(num_of_coordinates_in_vector):
            centroids[i][j] = "%.4f" % centroids[i][j]                                                                                   
        cent = ",".join(centroids[i])
        file.write(cent + "\n")

    # str_centroids = ",".join(centroids)
    # file.write(str_centroids)
    file.close()



main()