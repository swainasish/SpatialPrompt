"""
Created on Tue Jun 28 12:11:26 2022
@author: swainkasish@gmail.com
"""
#%% import libraries 
import numpy as np
import pandas as pd 
import time 
from sklearn.preprocessing import  StandardScaler,Normalizer,MinMaxScaler
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import KNeighborsRegressor
from alive_progress import alive_bar
from scipy.spatial import cKDTree
from  scipy.spatial.distance import cdist
from sklearn.cluster import KMeans,AgglomerativeClustering
from sklearn.linear_model import Ridge,RidgeCV
#%% fast spatial_prompt 
class SpatialDeconvolution:
    def __init__(self):
        pass
    
    def __preprocessing__(self,sc_array,st_array,sc_genes,st_genes,n_hvgs=1000):
        
           
        #Find common genes between SC and ST data 
        with alive_bar(6,title="Preprocessing Datasets :") as bar:
            bar()
            sc_genes,st_genes = np.array(sc_genes),np.array(st_genes)
            sc_df =  pd.DataFrame(sc_array,columns=sc_genes).astype("float32")
            st_df =  pd.DataFrame(st_array,columns=st_genes).astype("float32")
            bar()
            sc_df = sc_df.loc[:,~sc_df.columns.duplicated()]
            st_df = st_df.loc[:,~st_df.columns.duplicated()]
            bar()
            common_genes = np.intersect1d(sc_df.columns,st_df.columns)
            sc_cg_df = sc_df.loc[:,common_genes]
            n_cells,n_cg = sc_cg_df.shape
            if n_cg < n_hvgs:
                bar()
                sc_df = sc_df.loc[:,common_genes]
                st_df = sc_df.loc[:,common_genes]
                bar()
                bar()
            else: 
                if n_cells > 5000:
                    random_index = np.random.choice(n_cells,5000)
                    sc_cg_df = sc_cg_df.iloc[random_index,:]
                else:
                    pass
                sc_cg_df = self.__norm_hvg__(sc_cg_df)
                sc_vars = sc_cg_df.var().to_numpy()
                bar()
                sorted_variance_indexes = np.argsort(sc_vars)[::-1]
                bar()
                top_hvg_index = sorted_variance_indexes[0:n_hvgs]
                
                self.hvgs = sc_cg_df.columns[top_hvg_index]
                sc_df = sc_df.loc[:,self.hvgs]
                st_df = st_df.loc[:,self.hvgs]
                bar()
        n_hvg_final = sc_df.shape[1]
        print(f"{len(common_genes)} Common Genes Found {n_hvg_final} HVGs Retained")

        
        return sc_df,st_df
    
    def __random_spot_generator__(self,sc_arr_pro,sc_cell_type_labels,min_cell=5,max_cell=10):
        sc_arr_pro = np.array(sc_arr_pro)
        sc_cell_type_labels = np.array(sc_cell_type_labels)
        n_cell, n_gene = sc_arr_pro.shape
        
        #generate random cell number 
        random_cell_number = np.random.randint(low = min_cell, high=max_cell+1)
        
        #pick random rows in sc_data 
        
        random_indexes = np.random.choice(n_cell,random_cell_number)  
        random_sc_arrays = sc_arr_pro[random_indexes,:]
        random_sc_arrays = random_sc_arrays.sum(axis=0)
        
        random_sc_labels = sc_cell_type_labels[random_indexes]
        ct,ct_count = np.unique(random_sc_labels,return_counts=True)
        ct_perct = ct_count/np.sum(ct_count)
        label_dic = {}
        for i in range(len(ct)):
            label_dic[ct[i]]=ct_perct[i]
    
        
        return random_sc_arrays,label_dic
    
        
 
    
    def __simulate_st_spots__(self,sc_arr_pro,sc_cell_type_labels,min_cell=5,max_cell=10,
                              spot_ratio = [0.33,0.33,0.33]):
        
        """ This function generate spatial_spots from the single cell data 
            Using three scenarios mentioned in the paper 
            
            Parameters
            -------------------------------------------------------------------
            
            | sc_arr_pro : single data output from preprocessing function
            | sc_cell_type_labels : Cell type annotations of single cell data (now of sc_data rows = length of labels)
            | min_cell : (default = 5) Minium number of cell a spot have 
            | max_cell : (default = 10) Maxium number of cell a spot have 
            
            
        """

        # Convert datas to numpy array 
        sc_arr_pro = np.array(sc_arr_pro)
        sc_cell_type_labels = np.array(sc_cell_type_labels)
        n_cell, n_gene = sc_arr_pro.shape
        n_spot_to_simulated = 25586
        spot_ratio = np.array(spot_ratio)
        cri_1,cri_2,cri_3 = np.array((spot_ratio/sum(spot_ratio)) * n_spot_to_simulated,
                                     dtype=np.int64)
        n_spot_to_simulated = sum([cri_1,cri_2,cri_3])
        unique_cell_types = np.unique(sc_cell_type_labels)

        simulated_st_array = np.zeros([n_spot_to_simulated,n_gene])
        simulated_ct_prop = pd.DataFrame(np.zeros([n_spot_to_simulated,len(unique_cell_types)]))
        simulated_ct_prop.columns = unique_cell_types
        
        
        count=0
        with alive_bar(int(n_spot_to_simulated),title="Simulate Spatial Spots :") as bar:
            #Scenario.1 - Each spot having only one celltype
            for i in range(cri_1):
                random_cell_number = np.random.randint(low = min_cell, high=max_cell+1)
                #random_cell_type_pick
                random_cell_type = np.random.choice(unique_cell_types)
                cell_type_indexes = np.where(sc_cell_type_labels==random_cell_type)[0]
                
                #pick the indexes from the sc data
                index_to_pick = np.random.choice(cell_type_indexes,random_cell_number)
                random_array_frame = sc_arr_pro[index_to_pick,:]
                spot_1 = random_array_frame.sum(axis=0)
                simulated_st_array[count,:] = spot_1
                simulated_ct_prop.loc[count,random_cell_type] = 1  
                bar()
                count +=1
                
            # Scenario.2 - Each spot comprise of 2 cell type 
            if max_cell <= 2:   #designed for slideseq
                for i in range(cri_2):
                    sc_array , ct_dic = self.__random_spot_generator__(sc_arr_pro,sc_cell_type_labels,min_cell=0,max_cell=max_cell)
                    simulated_st_array[count,:] = sc_array
                    for i in ct_dic.keys():
                        simulated_ct_prop.loc[count,i]=ct_dic[i]
                    count=count+1
                    bar()
            else:
                for i in range(cri_2):
                    ct1,ct2 = np.random.choice(unique_cell_types,2)
                    ct1_array_index = np.where(ct1==sc_cell_type_labels)[0]
                    ct2_array_index = np.where(ct2==sc_cell_type_labels)[0]
                    
                    random_cell_number = np.random.randint(low = min_cell, high=max_cell+1)
                    ct1_cell_num = int(random_cell_number/2)
                    ct2_cell_num = int(random_cell_number - ct1_cell_num)
                    
                    ct1_index_to_pick = np.random.choice(ct1_array_index,ct1_cell_num)
                    ct2_index_to_pick = np.random.choice(ct2_array_index,ct2_cell_num)
                    concate_indexes = [*ct1_index_to_pick,*ct2_index_to_pick]
                    
                    array_frame = sc_arr_pro[concate_indexes,:]
                    array_frame = array_frame.sum(axis=0)
                    ct_prop = {}
                    ct_prop[ct1] = ct1_cell_num / random_cell_number
                    ct_prop[ct2] = ct2_cell_num / random_cell_number  
                    
                    simulated_st_array[count,:] = array_frame
                    for i in ct_prop.keys():
                        simulated_ct_prop.loc[count,i]=ct_prop[i]
                    
                    count= count+1
                    bar()
                    
            #Scenario-3 Each spot comprise of random cell type composition 
            for i in range(cri_3):
                random_cell_number = np.random.randint(low = min_cell, high=max_cell+1)
                random_indexes = np.random.choice(n_cell,random_cell_number)
                random_sc_df = sc_arr_pro[random_indexes,:]
                random_sc_label=sc_cell_type_labels[random_indexes]
                
                random_sc_df_sum = random_sc_df.sum(axis=0)
                ct,ct_count = np.unique(random_sc_label,return_counts=True)
                ct_count = ct_count/sum(ct_count)
                
                simulated_st_array[count,:]=random_sc_df_sum
                
                for i in range(len(ct)):
                    simulated_ct_prop.loc[count,ct[i]]=ct_count[i]
                
                count = count + 1
                bar()
                          
            
        print(f"{n_spot_to_simulated} Spatial Spots Simulated")
        
        return simulated_st_array,simulated_ct_prop
    
    
    def __normalisation1__(self,df):


        df=pd.DataFrame(df)
        
        df_row_sum = df.sum(axis=1)    #tpm
        df = df.div(df_row_sum,axis=0)* 10**6
        df = df.fillna(0)   #If a cell or spot having all counts are zero the infite values replaced by 0
        df = np.array(df)
        return df
    def __norm_hvg__(self,df):
        df_row_sum = df.sum(axis=1)    #tpm
        df = df.div(df_row_sum,axis=0)* 10**6
        df = df.fillna(0)   #If a cell or spot having all counts are zero the infite values replaced by 0
        return df
    
    
    def __normalisation2__(self,df):
        df=pd.DataFrame(df)
        
        df_row_sum = df.sum(axis=1)    #tpm
        df = df.div(df_row_sum,axis=0)* 10**6
        df = df.fillna(0)   #If a cell or spot having all counts are zero the infite values replaced by 0
        std = StandardScaler()
        norm = Normalizer()
        df = std.fit_transform(df)
        df = norm.fit_transform(df)
        return np.array(df)

    def __chunk_df__(self,df,label,n):
        n_cell,n_gene = df.shape
        df = np.array(df)
        label = np.array(label)
        if n_cell < n:
            pass
        else:
            random_points = np.random.randint(0,n_cell,n)
            df = df[random_points,:]
            label= label[random_points,:]
        return df,label
    def __domain_neighbor_index__(self,expr_matrix,x_cor,y_cor,n_neighbor=50):
        x_cor,y_cor = np.array(x_cor),np.array(y_cor)
        cor_mat = np.array([x_cor,y_cor],dtype=np.float32).T
        tree = cKDTree(cor_mat)    
        distances, nei_indexes = tree.query(cor_mat, k=n_neighbor)
        nei_indexes = nei_indexes[:,1:]
        n_spot,n_gene = expr_matrix.shape
        domain_index = np.zeros([n_spot,int(n_neighbor/2)])
        for i in range(n_spot):
            own_expression = expr_matrix[i,:]
            spot_nei_index = nei_indexes[i]
            nei_expression = expr_matrix[spot_nei_index,:]
            similarity  = 1-cdist(own_expression.reshape(1,-1),nei_expression,'cosine')[0]
            spot_dom_index = spot_nei_index[np.argsort(similarity)[::-1][0:int(n_neighbor/2)]]
            domain_index[i,:] = spot_dom_index
        domain_index = np.array(domain_index,dtype=np.int32)
        return domain_index
    
    def __domain_inspired__(self,expr_matrx,domain_index,n_itr=3,W=0.3):
        n_spot,n_gene = expr_matrx.shape
        expr_spatial_neighbor = np.array([expr_matrx[domain_index[i,:],:].mean(axis=0) for i in range(n_spot)],dtype=np.float32)
        for weight in range(n_itr):
            nei_weightage = W/(weight+1)
            old_weightage = 1-nei_weightage
            expr_micro_env = np.array([expr_spatial_neighbor[domain_index[i,:],:].mean(axis=0) for i in range(n_spot)],dtype=np.float32)
            expr_spatial_neighbor = (old_weightage * expr_spatial_neighbor) + (nei_weightage * expr_micro_env)
        return expr_spatial_neighbor
    
    
    def __spatial_integration__(self,simu_arr,real_arr,x_cor,y_cor,n_neighbor =40,n_itr=3):
        with alive_bar(5,title="Capturing spatial microenvironment relation:") as bar:
            simu_arr_norm,real_arr_norm = self.__normalisation1__(simu_arr),self.__normalisation1__(real_arr)
            n_spot,n_gene = real_arr_norm.shape
            domain_index = self.__domain_neighbor_index__(real_arr_norm, x_cor, y_cor,n_neighbor =n_neighbor)
            std_scale = StandardScaler()
            bar()
            real_arr_norm_std = std_scale.fit_transform(real_arr_norm)
            simu_arr_norm_std = std_scale.fit_transform(simu_arr_norm)
            bar()
            real_exp_spatial = self.__domain_inspired__(real_arr_norm, domain_index,n_itr=n_itr )
            ridge_cv = RidgeCV(alphas=[1,10,25,50,100,500,1000,5000])
            bar()
            chunk_st, chunk_sc = self.__chunk_df__(real_arr_norm_std,real_exp_spatial,5000)
            ridge_cv.fit(chunk_st, chunk_sc)
            ridge_model = Ridge(alpha=ridge_cv.alpha_)
            bar()
            ridge_model.fit(real_arr_norm_std,real_exp_spatial)
            simu_exp_spatial = ridge_model.predict(simu_arr_norm_std)
            simu_exp_spatial[simu_exp_spatial<0]=0            
            real_nei_arr = np.hstack([real_arr_norm,real_exp_spatial])
            simu_nei_arr = np.hstack([simu_arr_norm,simu_exp_spatial])
            bar()   
        return real_nei_arr,simu_nei_arr

    def __fit_predict__(self,real_st,simu_st,simu_prop,return_prop = True):
        
        """
        Predict the cell type proportations from the Spatial data 
        
        Parameters
        ----------------------------------------------------------------------
        st_data : st_data output from the preprocessing function 
        return_low_exp_celltypes  : low expressed cell proportions needed or not
        cell_prop_threshold : Threshold for low expression celltype
            
        """
        with alive_bar(5,title="Spot Denvolution:") as bar:
            simu_st_norm = self.__normalisation2__(simu_st)
            real_st_norm = self.__normalisation2__(real_st)
            bar()
      
            model = KNeighborsRegressor(n_neighbors=250)
            bar()
            # print(best_alpha)
            model.fit(simu_st_norm,simu_prop)
            bar()
            predict = model.predict(np.array(real_st_norm))
            bar()
            predict[predict<0.0] = 0 
            row_sum = predict.sum(axis=1).reshape(-1,1)
            predict = predict / row_sum
            predict = pd.DataFrame(predict,
                                    columns=simu_prop.columns)
            bar()
            if return_prop==False:
                predict = pd.Categorical(predict.idxmax(axis=1))

        return predict
    
    
    def predict_cell_prop(self,sc_array,st_array,
                          sc_genes,st_genes,sc_labels,
                          x_cord,y_cord,
                          n_hvgs=1000,min_cell=10,max_cell=15,
                          return_prop=True,spot_ratio= [0.33,0.33,0.33],n_neighbor=45,n_itr=3):
        '''
        ### Description
        This program perform spot deconvolution in spatial data using scRNA-seq data reference.
        
        ### Parameters
        - `sc_array`: Matrix of Single cell data, where rows are the cells and columns are the genes.
        - `st_array`: Matrix of Spatial data, where rows are the cells and columns are the genes.
        - `sc_genes`: Gene names of the `sc_array` matrix.
        - `st_genes`: Gene names of the `st_array` matrix.
        - `sc_labels`: Cell type proportions of `sc_array`.
        - `x_cord`: X coordinate array of spatial data.
        - `y_cord`: Y coordinate array of spatial data.
        - `n_hvgs` (default=1000): Number of high variance genes to consider for analysis.
        - `min_cell` (default=10): Minimum number of cells to simulate the spatial spot.
        - `max_cell` (default=15): Maximum number of cells to simulate the spatial spot.
        - `return_prop` (default=True): Return proportions of cell types if true, else return the cell type having a higher proportion.
        - `spot_ratio` (default=[0.33, 0.33, 0.33]): Ratio of proportions of spots to be simulated using criteria 1/2/3 mentioned in the paper. If the labels are ambiguous cell types (e.g., EX_L3_4_5 have cell types of L3 AND L4), then `spot_ratio` should be provided as a list, e.g., `[0, 0, 1]`.
        - `n_neighbor` (default=45): Number of neighbors to consider for weighted mean expression calculation.
        - `n_itr` (default=3): Number of iterations message passing layer pull information from neighbors.
        '''
        t1 = time.time()
        sc_processed,st_processed = self.__preprocessing__(sc_array, st_array,
                                                           sc_genes, st_genes,n_hvgs=n_hvgs) 

        simu_st_arr,simu_ct_prop = self.__simulate_st_spots__(sc_processed,sc_labels,min_cell=min_cell,
                                                               max_cell=max_cell,spot_ratio=spot_ratio)
        
        st_processed,simu_st_arr = self.__spatial_integration__(simu_st_arr,st_processed, x_cord, y_cord,n_neighbor=n_neighbor,
                                                                n_itr=n_itr)

        st_pred_ct = self.__fit_predict__(st_processed,simu_st_arr,simu_ct_prop,return_prop=return_prop)

        t2= time.time()
        print(f"Total Time spent: {t2-t1} Sec")
    
        return st_pred_ct
    def cluster_cell_annotation(self,ct_deconv,cluster_annot,threshold=0.30):
        dicti_annot = {}
        cluster_nums = len(np.unique(cluster_annot))
        for i in range(cluster_nums):
            clus_temp1_index = cluster_annot == i
            clus_deconvo = ct_deconv[clus_temp1_index].mean(axis=0)
            clus_deconvo_major = clus_deconvo[clus_deconvo>threshold]
            if len(clus_deconvo_major)==0:
                thres = threshold - 0.10
                clus_deconvo_major = clus_deconvo[clus_deconvo>thres]
            dicti_annot[i] = ' / '.join(list(clus_deconvo_major.index))
            if len(clus_deconvo_major)==0:
                dicti_annot[i] = "ambiguous"
        mod_cluster_annot = pd.Categorical([dicti_annot[i]for i in cluster_annot])
        return dicti_annot,mod_cluster_annot
    

#%% 
class SpatialCluster:
    def __init__(self):
        pass
    def __hvg_detect__(self,df_array,N=1000):
        df_array = np.array(df_array)
        st_vars = df_array.var(axis=0)
        top_n_index = np.argsort(st_vars)[::-1][0:N]
        df_reduced = df_array[:,top_n_index] 
        return df_reduced
    
    def __kmean_clus__(self,df,n_clus="auto"):
        # std = StandardScaler()
        std = MinMaxScaler()
        pca =   TruncatedSVD(n_components=50)
        st_df_graph_std = std.fit_transform(df)
        st_df_graph_pca = pca.fit_transform(st_df_graph_std)
        if n_clus=="auto":
            n_clus = self.__auto_k_find__(st_df_graph_pca)
        kmean = KMeans(n_clusters=n_clus,n_init=20)
        kmean.fit(st_df_graph_pca)
        labels = np.array(kmean.labels_)
        return labels
    
    def __auto_k_find__(self,data):
        n_row,_ = data.shape
        random_row = np.random.choice(n_row,5000)
        data_chunk = data[random_row,:]
        agglo_clust = AgglomerativeClustering(n_clusters=None,distance_threshold=15)
        agglo_clust.fit(data_chunk)
        labels = agglo_clust.labels_
        best_k = len(np.unique(labels))
        
        return best_k
    
    def __normalisation3__(self,df):
        df=pd.DataFrame(df)     
        df_row_sum = df.sum(axis=1)    #tpm
        df = df.div(df_row_sum,axis=0)* 10**6
        df = df.fillna(0)   
        minmax= MinMaxScaler()
        df = minmax.fit_transform(df)
        return np.array(df)
    
    def fit_predict(self,st_array,x_cord,y_cord,n_neighbor=20,n_itr=3,n_cluster="auto",W=0.4,n_hvgs=1000):
        """
        ### Description
        This program perform spatial clustering for spatial data .
        
        ### Parameters
        - `st_array`: Matrix of Spatial data, where rows are the cells and columns are the genes.
        - `x_cord`: X coordinate array of spatial data.
        - `y_cord`: Y coordinate array of spatial data.
        - `n_hvgs` (default=1000): Number of high variance genes to consider for analysis.
        - `n_neighbor` (default=45): Number of neighbors to consider for weighted mean expression calculation.
        - `n_itr` (default=3): Number of iterations message passing layer pull information from neighbors.
        """
        t1 = time.time()
        with alive_bar(5,title="Spatial clustering: ") as bar:
            st_array_hvg = self.__hvg_detect__(st_array,N=n_hvgs)
            bar()
            st_array_hvg_norm = self.__normalisation3__(st_array_hvg)
            deconv_class = SpatialDeconvolution()
            domain_index = deconv_class.__domain_neighbor_index__(st_array_hvg_norm, x_cord, y_cord,n_neighbor=n_neighbor)
            bar()
            st_array_domain= deconv_class.__domain_inspired__(st_array_hvg_norm, domain_index,n_itr=n_itr,W=W )
            bar()
            st_array_hvg_norm_domain = np.hstack([st_array_hvg_norm,st_array_domain])
            bar()
            cluster_assignment = self.__kmean_clus__(st_array_hvg_norm_domain,n_clus = n_cluster)
            cluster_assignment = pd.Categorical(cluster_assignment)
            bar()
        t2 = time.time()
        print(f"Executed in {t2-t1} second")
        
        return cluster_assignment
#%% SpatialPromptDB
class Database_linker:
    def __init__(self):
        pass
    def available_tissues(self,species="human"):
        pass
    def reference_download(species="human",tissue="cortex"):
        pass
        





