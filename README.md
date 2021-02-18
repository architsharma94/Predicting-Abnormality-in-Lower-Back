# Predicting-Abnormality-in-Lower-Back

## Repository guide:
* Train Data => train_spine.csv
* Test Data => test_spine.csv
* Analysis and Reporting => Predicting Abnormality in Lower Back.docx
* Code => code.R
* Test runs => run-1 to 5.RData

## Introduction

With the onset of the COVID-19 pandemic, work and studies have shifted online with working on laptops from home. Being in lockdown and studies being remote students have started spending long hours in front of laptops camped out on bed, kitchen countertops or tables. Due to bad postures, people are experiencing new aches and physical strain especially lower back pains more commonly.  This ignited our curiosity on whether people with certain body structure are more susceptible to have lower back pain and led to the topic for this project.
The bottom part of the back or lumbar spine has network of spinal muscles, interconnected bones, nerves and is connected to the pelvis, working in tandem. Generating lots of movements and bearing the weight of the upper body, the lower back becomes vulnerable to pain. Though one of the most common type of pain, lumbago (lower back pain) can be acute or chronic. It could be indications of degenerative or nerve and spinal problems. (Akhil Chhatre, n.d.) (John Peloza, 2017)
Lower back pain may also be caused by damages in joints, ligaments or degenerating intervertebral disc or irritation in the different components of the back. In this project, we differentiate body structures as being normal and abnormal, thereby determining susceptibility to lower back pain, based on the physical attributes of the spine and pelvis. With the information contained in the independent features or predictor variables, and the dependent variable being binary, the focus of this project would be to classify or predict the type of body structure of person using robust Bayesian binary Logistic Regression analysis methods. (Stroke, n.d.)

## Conclusion

To conclude, a robust logistic regression model for classifying the condition of an individual’s lower back was successfully built. Feature selection was performed to reduce the dimensions of model by using principles of correlation and statistical insignificance. 
As per the results obtained, the model simulated in 6th run was considered ideal for prediction. The selection was based on principle of parsimony and efficiency. Of the 3 models considered for gauging predictive power, the configuration from 6th run took lesser computational time than 7th run and utilized lesser number of predictors than 5th configuration.
The reduction of dimensions from 5th run to 6th run, yielded similar results in terms of predictive power. Hence, it was also inferred that the MCMC module automatically disregards the insignificant parameters for prediction.
The final predictive model is given below:

Y = -14.2 - 0.079*Pelvic_Tilt + 0.093*Sacral_Slope+ 0.107*Pelvic_Radius - 0.147*Degree_Spondylolisthesis

For a unit increase in pelvic tilt and degree of spondylolisthesis, the log odds of outcome would decrease by 0.079 and 0.147, respectively. On the contrary, for a unit increase in sacral slope and pelvic radius, the log odds of outcome would increase by 0.093 and 0.107, respectively.
It can be observed for 1st prediction that the probability of normal lower back was 0.00005 (~0). This probability extends from 0 to 0.00099. This means that the probability of lower back being abnormal is very high. Similar inferences can be drawn for predictions 1 to 7. However, in case of 8th prediction, the probability of normal lower back was 0.65, which is well above the 0.5 threshold. This means that the individual had a higher probability of having a normal lower back. 
The model was able to predict the condition of lower back of an individual with 90% accuracy. 

