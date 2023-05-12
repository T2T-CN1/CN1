from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
import sys

def load_data(file_path):
    '''
    readid, genotype
    '''
    with open(file_path,'r') as file:
        lines = file.readlines()
    read_name = []
    nucCodes = []
    for line in lines:
        items = line.strip().split('\t')
        if line.startswith('Reads'):continue
        read_name.append(items[0])
        nucCodes.append([items[i] for i in range(1,len(items))])

    return read_name,nucCodes


def cluster(data,read_name,cluster_num):
    '''
    构造KMeans聚类，返回聚类簇和对应的族中心坐标加和
    '''
    km = KMeans(n_clusters = cluster_num)
    y_kmeans = km.fit_predict(data)
    # predict cluster number
    wcss = []
    for i in range(1,11):
        kmeans = KMeans(n_clusters = i, max_iter = 300, n_init = 10, init = 'k-means++', random_state = 0)
        kmeans.fit(data)
        wcss.append(kmeans.inertia_)
    #plt.plot(range(1,11), wcss)
    #plt.title('The Elbow Method')
    #plt.xlabel('Number of Clusters')
    #plt.ylabel('WCSS')
    #plt.show()

    # Visualizing the clusters
    fig, axs = plt.subplots(figsize=(15, 5))
    centers = kmeans.cluster_centers_
    #ic(data[:,0])
    #axs.scatter(data[:, 0], data[:, 1], s=10, c=kmeans.labels_)
    #axs.scatter(centers[:, 0], centers[:, 1], c="r", s=20)
    #plt.scatter(centers[:, 0],  kmeans.cluster_centers_[:, 1], c = 'yellow', label = 'Centroids')

    #plt.show()


    # 每个簇中心的坐标 cluster_centers_,
    rdna_profile = np.sum(km.cluster_centers_,axis=1)

    # 构造存放cluster_num个聚类的列表
    read_cluster = []
    for i in range(cluster_num):
        read_cluster.append([])

    for i in range(len(read_name)):
        read_cluster[y_kmeans[i]].append(read_name[i])

    return read_cluster,rdna_profile



def main():
    file_path = sys.argv[1]
    cluster_num = 5
    read_name,nucCodes = load_data(file_path)
    read_cluster,rdna_profile = cluster(nucCodes,read_name,cluster_num)

    for i in range(len(read_cluster)):
        #print("distance:%.2f" %rdna_profile[i])
        #print(read_cluster[i])
        index = i + 1
        for read in read_cluster[i]:
            print(f"cluster{index}\t{rdna_profile[i]}\t{read}")

if __name__== '__main__':
    if len(sys.argv) < 2:
        sys.exit(f"python3 {sys.argv[0]} *.matrix")
    else:
        main()
