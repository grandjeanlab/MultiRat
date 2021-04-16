MultiRAT environement preparation
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

# Foreword

This and the follow are jupyter files which contains all the code for
reproducing my analysis and detail the process. The code is meant to be
followed step-wise. The raw fMRI dataset will not be publicly available
before the project preprint publication on BioRxiv. The raw fMRI dataset
can be made available prior to publication upon request and review from
the authors.

If re-using some of the scripts, please follow citations guidelines for
the software used. I’ve provided the links to the software wherever
possible. See also the [license](../LICENSE.md) for this software.

The code is executed in `bash` (RABIES preprocessing) and `python` (analysis
and plots).

The `python` libraries are detailed in a [yml file](../environment.yml).

To reproduce the code contained within this software, please follow
these steps  
1\. get the required dependencies 
2\. Install python packages (`conda env create -f environment.yml`)
3\. Build RABIES singularity image (`singularity build rabies-0.2.1.simg docker://gabdesgreg/rabies:0.2.1`)
4\. Update the variables `init_folder` and `analysis_folder` at the top of the notebooks. 



```python
# Python environement
! python --version
! conda --version
! jupyter --version
```

    Python 2.7.5
    conda 4.9.0
    jupyter core     : 4.6.3
    jupyter-notebook : 6.0.3
    qtconsole        : 4.7.5
    ipython          : 7.16.1
    ipykernel        : 5.3.2
    jupyter client   : 6.1.6
    jupyter lab      : 2.1.5
    nbconvert        : 5.6.1
    ipywidgets       : 7.5.1
    nbformat         : 5.0.7
    traitlets        : 4.3.3



```python
# FSL version 
! flirt -version
```

    FLIRT version 6.0



```python
# install python packages
! pip install nipype pandas numpy seaborn matplotlib nibabel 
```

    Defaulting to user installation because normal site-packages is not writeable
    Collecting nipype
      Downloading nipype-1.6.0-py3-none-any.whl (3.1 MB)
    [K     |████████████████████████████████| 3.1 MB 9.1 MB/s eta 0:00:01
    [?25hRequirement already satisfied: pandas in /opt/anaconda3/2020.07/lib/python3.8/site-packages (1.0.5)
    Requirement already satisfied: numpy in /opt/anaconda3/2020.07/lib/python3.8/site-packages (1.18.5)
    Requirement already satisfied: seaborn in /opt/anaconda3/2020.07/lib/python3.8/site-packages (0.10.1)
    Requirement already satisfied: matplotlib in /opt/anaconda3/2020.07/lib/python3.8/site-packages (3.2.2)
    Collecting nibabel
      Using cached nibabel-3.2.1-py3-none-any.whl (3.3 MB)
    Requirement already satisfied: networkx>=2.0 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from nipype) (2.4)
    Collecting etelemetry>=0.2.0
      Downloading etelemetry-0.2.2-py3-none-any.whl (6.2 kB)
    Requirement already satisfied: packaging in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from nipype) (20.4)
    Requirement already satisfied: python-dateutil>=2.2 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from nipype) (2.8.1)
    Collecting simplejson>=3.8.0
      Downloading simplejson-3.17.2-cp38-cp38-manylinux2010_x86_64.whl (137 kB)
    [K     |████████████████████████████████| 137 kB 52.9 MB/s eta 0:00:01
    [?25hCollecting pydot>=1.2.3
      Downloading pydot-1.4.2-py2.py3-none-any.whl (21 kB)
    Collecting rdflib>=5.0.0
      Using cached rdflib-5.0.0-py3-none-any.whl (231 kB)
    Requirement already satisfied: filelock>=3.0.0 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from nipype) (3.0.12)
    Requirement already satisfied: scipy>=0.14 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from nipype) (1.5.0)
    Requirement already satisfied: click>=6.6.0 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from nipype) (7.1.2)
    Collecting prov>=1.5.2
      Downloading prov-2.0.0-py3-none-any.whl (421 kB)
    [K     |████████████████████████████████| 421 kB 58.2 MB/s eta 0:00:01
    [?25hCollecting traits!=5.0,>=4.6
      Downloading traits-6.2.0-cp38-cp38-manylinux2010_x86_64.whl (5.1 MB)
    [K     |████████████████████████████████| 5.1 MB 45.8 MB/s eta 0:00:01     |█████████████████████████████▏  | 4.6 MB 45.8 MB/s eta 0:00:01
    [?25hRequirement already satisfied: pytz>=2017.2 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from pandas) (2020.1)
    Requirement already satisfied: kiwisolver>=1.0.1 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from matplotlib) (1.2.0)
    Requirement already satisfied: cycler>=0.10 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from matplotlib) (0.10.0)
    Requirement already satisfied: pyparsing!=2.0.4,!=2.1.2,!=2.1.6,>=2.0.1 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from matplotlib) (2.4.7)
    Requirement already satisfied: decorator>=4.3.0 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from networkx>=2.0->nipype) (4.4.2)
    Requirement already satisfied: requests in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from etelemetry>=0.2.0->nipype) (2.24.0)
    Collecting ci-info>=0.2
      Downloading ci_info-0.2.0-py3-none-any.whl (6.9 kB)
    Requirement already satisfied: six in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from packaging->nipype) (1.15.0)
    Collecting isodate
      Using cached isodate-0.6.0-py2.py3-none-any.whl (45 kB)
    Requirement already satisfied: lxml>=3.3.5 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from prov>=1.5.2->nipype) (4.5.2)
    Requirement already satisfied: chardet<4,>=3.0.2 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from requests->etelemetry>=0.2.0->nipype) (3.0.4)
    Requirement already satisfied: idna<3,>=2.5 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from requests->etelemetry>=0.2.0->nipype) (2.10)
    Requirement already satisfied: certifi>=2017.4.17 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from requests->etelemetry>=0.2.0->nipype) (2020.6.20)
    Requirement already satisfied: urllib3!=1.25.0,!=1.25.1,<1.26,>=1.21.1 in /opt/anaconda3/2020.07/lib/python3.8/site-packages (from requests->etelemetry>=0.2.0->nipype) (1.25.9)
    Installing collected packages: ci-info, etelemetry, simplejson, pydot, isodate, rdflib, prov, nibabel, traits, nipype
    Successfully installed ci-info-0.2.0 etelemetry-0.2.2 isodate-0.6.0 nibabel-3.2.1 nipype-1.6.0 prov-2.0.0 pydot-1.4.2 rdflib-5.0.0 simplejson-3.17.2 traits-6.2.0



```python
! conda env export > ../environment.yml
```
