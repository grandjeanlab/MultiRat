MultiRAT analysis code
================
Joanes Grandjean

![rat art](../assets/img/rat_art.png)

In this section, I provide a brief summary of the datasets that were
supplied for this study.

Normalization steps: - Columns with .orig have been left unaltered from
how it was provided by the dataset owner. - Strain was simplified. - Age
normalized to months, in 2 month bins. When long range are provided, the
middle bin is selected. - Weight normalized to in 50 g bins. When long
range are provided, the middle bin is selected. - Dexmedetomidine dose
converted to medetomidine (\*0.5) for simplicity. - Injectable
anesthesia dose converted to mg/kg - Bolus injectable added to induction
column, infusions kept in maintenance column. - Multiple anesthesia
agents for either induction or maintenance indicated with ‘/’.
Corresponding doses follow same order. - Measurement post induction,
breathing rate, heart rate rounded to the nearest tenth - TR and TE
converted to s - Bandwidth converted to Hz

First, I examine subject distribution

``` r
library(ggplot2)
df <- read.csv2('../assets/table/meta_data.tsv',sep='\t')
df$rat.age <- factor(df$rat.age, levels = c("0-2", "2-4", "4-6","6-8","8-10","10-12","12-14","14-16","16-18","18-20"))
```

``` r
summary(df$rat.strain)
```

    ##    Fischer 344  Lister Hooded     Long Evans Sprague Dawley         Wistar 
    ##             70             20             70            145            136

Our collection of datasets is enriched in Wistar and Sprage Dawley rats,
with Lister Hooded being least frequently used.

``` r
summary(df$rat.sex)
```

    ## Female   Male 
    ##    130    311

Consistent with the “male bias” in science, we nearly have 1:3
Female/Male ratio.

``` r
#summary(df$rat.age)
ggplot(df, aes(x = rat.age, fill=rat.strain)) + geom_dotplot(binwidth = 0.1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](proj_dataset_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> We tend
to have younger rats, aged 0-6 months. Older rats are also available,
but strain distribution is not even.

``` r
ggplot(df, aes(x = rat.weight, fill=rat.strain)) + geom_dotplot(binwidth = 0.15) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](proj_dataset_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> Weight
distribution as a function of strain. Oddly, Fisher, despite being the
oldest, are the lightest. The rats are mostly in the 300-350 range,
which is consistent with my expectations.

``` r
ggplot(df, aes(x = as.factor(anesthesia.breathing.rate), fill=rat.strain)) + geom_dotplot(binwidth = 0.15) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](proj_dataset_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Selecting first rat of each dataset to get a representation of dataset parameter representations. 

df.sub <- df[df$rat.ses == 1 & sprintf('%02d', df$rat.sub %% 100) == "00", ]
summary(df.sub[,c(20,21,23,28,36,38,42,46, 47, 48, 49)])
```

    ##             exp.type                     anesthesia.induction
    ##  resting-state  :28   awake                        : 1       
    ##  stimulus-evoked:10   isoflurane                   :19       
    ##                       isoflurane / alpha-chloralose: 1       
    ##                       isoflurane / medetomidine    :11       
    ##                       ketamine / xylazine          : 2       
    ##                       urethane                     : 4       
    ##                                                              
    ##                anesthesia.maintenance anesthesia.breathing.assistance
    ##  alpha-chloralose         : 2         free-breathing:32              
    ##  awake                    : 1         ventilated    : 5              
    ##  isoflurane               :14         NA's          : 1              
    ##  isoflurane / medetomidine: 8                                        
    ##  medetomidine             : 8                                        
    ##  urethane                 : 5                                        
    ##                                                                      
    ##              MRI.vendor MRI.field.strength   anat.sequence func.sequence
    ##  Bruker           :30   11.1: 1            3D-FLASH : 1    GE-EPI:32    
    ##  Magnex Scientific: 1   14.1: 4            3D-RARE  : 1    SE-EPI: 6    
    ##  Varian/Agilent   : 7   4.7 : 3            bSSFP    : 1                 
    ##                         7   :14            FSE      : 3                 
    ##                         9.4 :16            MP2RAGE  : 1                 
    ##                                            RARE     :16                 
    ##                                            turboRARE:15                 
    ##     func.TR      func.TE     func.FA     
    ##  1      :13   0.015  :9   Min.   :15.00  
    ##  2      :12   0.018  :9   1st Qu.:60.00  
    ##  1.5    : 7   0.02   :4   Median :72.50  
    ##  1.6    : 2   0.045  :4   Mean   :73.88  
    ##  0.5    : 1   0.012  :3   3rd Qu.:90.00  
    ##  0.7    : 1   0      :2   Max.   :90.00  
    ##  (Other): 2   (Other):7   NA's   :4

Finally, I look at the dataset parameter distribution relevant for our
analysis. I obtained more resting-state than stimulus-evoked dataset
(3:1). Isoflurane is the most common way for induction, no surprise
here. Isoflurane, medetomidine, and their combination are also the most
common methods to maintain the animals in the scanner. No surprise.
Free-breathing, Bruker, 7 - 9.4T fields, and GE-EPI are the norms.
Again, no surprise. TR are in the 1-2 s range.

Overall, the dataset parameter distribution is aligned with our previous
observations about [fMRI in
rodents](https://www.frontiersin.org/articles/10.3389/fninf.2019.00078/full)

``` r
df.na <-df[!is.na(df$anesthesia.breathing.rate) & !is.na(df$rat.strain) & !is.na(df$rat.sex) & !is.na(df$rat.age) & !is.na(df$anesthesia.maintenance),]

mod.breathing.full <- lm(anesthesia.breathing.rate ~ rat.strain + rat.sex  + rat.age + anesthesia.maintenance, df.na,na.action = na.omit)

mod.breathing.strain <- update(mod.breathing.full, . ~ . - rat.strain)
anova(mod.breathing.full,mod.breathing.strain)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.breathing.rate ~ rat.strain + rat.sex + rat.age + 
    ##     anesthesia.maintenance
    ## Model 2: anesthesia.breathing.rate ~ rat.sex + rat.age + anesthesia.maintenance
    ##   Res.Df   RSS Df Sum of Sq      F  Pr(>F)  
    ## 1    311 40799                              
    ## 2    314 41902 -3   -1103.5 2.8039 0.03996 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod.breathing.sex <- update(mod.breathing.full, . ~ . - rat.sex)
anova(mod.breathing.full,mod.breathing.sex)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.breathing.rate ~ rat.strain + rat.sex + rat.age + 
    ##     anesthesia.maintenance
    ## Model 2: anesthesia.breathing.rate ~ rat.strain + rat.age + anesthesia.maintenance
    ##   Res.Df   RSS Df Sum of Sq      F  Pr(>F)  
    ## 1    311 40799                              
    ## 2    312 41497 -1    -698.4 5.3237 0.02169 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod.breathing.age <- update(mod.breathing.full, . ~ . - rat.age)
anova(mod.breathing.full,mod.breathing.age)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.breathing.rate ~ rat.strain + rat.sex + rat.age + 
    ##     anesthesia.maintenance
    ## Model 2: anesthesia.breathing.rate ~ rat.strain + rat.sex + anesthesia.maintenance
    ##   Res.Df   RSS Df Sum of Sq      F    Pr(>F)    
    ## 1    311 40799                                  
    ## 2    317 46802 -6   -6003.5 7.6273 1.185e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod.breathing.maintenance <- update(mod.breathing.full, . ~ . - anesthesia.maintenance)
anova(mod.breathing.full,mod.breathing.maintenance)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.breathing.rate ~ rat.strain + rat.sex + rat.age + 
    ##     anesthesia.maintenance
    ## Model 2: anesthesia.breathing.rate ~ rat.strain + rat.sex + rat.age
    ##   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
    ## 1    311  40799                                  
    ## 2    315 115672 -4    -74874 142.69 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

In an exploratory analysis, I test the effects associated with breathing
rates in rats. This analysis shows small effects of strain, sex, mild
effect of age, and very strong effect of anesthesia used for
maintenance. Both Urethane and medetomidine groups have higher breathing
rates.

``` r
ggplot(df.na, aes(y = as.factor(anesthesia.breathing.rate), x=anesthesia.maintenance, color=rat.strain)) + geom_jitter(size = 0.15) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](proj_dataset_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
df.na <-df[!is.na(df$anesthesia.heart.rate) & !is.na(df$rat.strain) & !is.na(df$rat.sex) & !is.na(df$rat.age) & !is.na(df$anesthesia.maintenance),]

mod.heart.full <- lm(anesthesia.heart.rate ~ rat.strain + rat.sex  + rat.age + anesthesia.maintenance, df.na,na.action = na.omit)


#summary(mod.breathing.full)

mod.heart.strain <- update(mod.heart.full, . ~ . - rat.strain)
anova(mod.heart.full,mod.heart.strain)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.heart.rate ~ rat.strain + rat.sex + rat.age + anesthesia.maintenance
    ## Model 2: anesthesia.heart.rate ~ rat.sex + rat.age + anesthesia.maintenance
    ##   Res.Df     RSS Df Sum of Sq      F    Pr(>F)    
    ## 1    158  103473                                  
    ## 2    160 1127666 -2  -1024192 781.95 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod.heart.sex <- update(mod.heart.full, . ~ . - rat.sex)
anova(mod.heart.full,mod.heart.sex)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.heart.rate ~ rat.strain + rat.sex + rat.age + anesthesia.maintenance
    ## Model 2: anesthesia.heart.rate ~ rat.strain + rat.age + anesthesia.maintenance
    ##   Res.Df    RSS Df Sum of Sq      F Pr(>F)
    ## 1    158 103473                           
    ## 2    159 104062 -1   -588.27 0.8983 0.3447

``` r
mod.heart.age <- update(mod.heart.full, . ~ . - rat.age)
anova(mod.heart.full,mod.heart.age)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.heart.rate ~ rat.strain + rat.sex + rat.age + anesthesia.maintenance
    ## Model 2: anesthesia.heart.rate ~ rat.strain + rat.sex + anesthesia.maintenance
    ##   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
    ## 1    158 103473                                  
    ## 2    161 381366 -3   -277893 141.44 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod.heart.maintenance <- update(mod.heart.full, . ~ . - anesthesia.maintenance)
anova(mod.heart.full,mod.heart.maintenance)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: anesthesia.heart.rate ~ rat.strain + rat.sex + rat.age + anesthesia.maintenance
    ## Model 2: anesthesia.heart.rate ~ rat.strain + rat.sex + rat.age
    ##   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
    ## 1    158 103473                                  
    ## 2    159 234121 -1   -130647 199.49 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

For heart rate, fewer data points are available. We do find the expected
medetomidine effect on heart rate, leading to overall lower heart rate
(Bradycardia). Because not all factors are orthogonal (e.g. strain and
anesthesia maintenance), effect interpretations are rendered difficult.
Similarly, effects of anesthesia dose, route, time, and induction are
likely beyond our reach.

``` r
ggplot(df.na, aes(y = as.factor(anesthesia.heart.rate), x=anesthesia.maintenance, color=rat.strain)) + geom_jitter(size = 0.15) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](proj_dataset_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
