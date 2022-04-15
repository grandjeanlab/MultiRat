import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.image as mpimg
import os
import re
import random
import glob

def close_event():
    plt.close() #timer calls this function after 3 seconds and closes the window 



ica_dir='/project/4180000.19/multiRat/export/aromas_qa/'
ica_files = glob.glob(ica_dir+'*')  
out_dir='/project/4180000.19/multiRat/export/aromas_qa_test/'
os.makedirs(out_dir, exist_ok=True)

ds_id=range(1001,1053)

for i in ds_id:
    r = re.compile(str(i))
    ica_select=list(filter(r.findall, ica_files))
    ica_select=random.sample(ica_select,2)
    img_cat = [''] * 10
    for j in ica_select:
        for k in range(1,10):
            fig = plt.figure(figsize=(6,12))
            timer = fig.canvas.new_timer(interval = 2000) #creating a timer object and setting an interval of 3000 milliseconds
            timer.add_callback(close_event)
            
            img = mpimg.imread(j+'/IC_'+str(k)+'_thresh.png')
            imgplot = plt.imshow(img)
            timer.start()
            plt.show()
            img_cat[k-1] = input("(s)ignal, (n)oise, (u)nsure \n")
    
        signal_cat = [x for x, z in enumerate(img_cat) if z == 's'] 
        noise_cat = [x for x, z in enumerate(img_cat) if z == 'n']
        unsure_cat = [x for x, z in enumerate(img_cat) if z == 'u']

        signal_cat = [x+1 for x in signal_cat]
        noise_cat = [x+1 for x in noise_cat]
        unsure_cat = [x+1 for x in unsure_cat]

        my_file = open(j+'/classification.txt', 'r')
        aroma_cat = my_file.read().split(',')
        aroma_cat = [int(x) for x in aroma_cat]

        false_neg=len(list(set(aroma_cat) & set(signal_cat)))
        true_neg=len(list(set(aroma_cat) & set(noise_cat)))
        unsure=len(list(set(aroma_cat) & set(unsure_cat)))

        with open(os.path.join(out_dir,os.path.basename(j)+'_false_neg.txt'), 'w') as f:
            f.write(str(false_neg))
        with open(os.path.join(out_dir,os.path.basename(j)+'_true_neg.txt'), 'w') as f:
            f.write(str(true_neg))
        with open(os.path.join(out_dir,os.path.basename(j)+'_unsure_neg.txt'), 'w') as f:
            f.write(str(unsure))

