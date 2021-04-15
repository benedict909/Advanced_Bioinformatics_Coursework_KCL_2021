---
title: "Advanced Bioinformatics 2021 assessment"
author: 'Candidate Number: 18040'
date: "23/03/2021"
output:
  word_document: default
  html_document: default
  pdf_document: default
---


## Task 3.1

Using the `sum()` function to evaluate the sum of all intergers between 5 and 55.
```{r sum}
sum(c(5:55))
```
<br>

## Task 3.2 

A function that calculates the sum of all integers between 5 and n. 

```{r}
sumfun = function(n){
  output = sum(c(5:n))
  return(output)
}
sumfun(10)
sumfun(20)
sumfun(100)
```
<br>

## Task 3.3 

Loop to calculate and print the first 5 terms of the Fibonacci sequence. 

```{r}
for(i in 1:12){
  if(i == 1){
    current_num = i
    prev_num = 0 # as no previous num in the sequence 
  }else{
    current_num = prev_num + prev_prev_num 
  }
    print(current_num) # print current entry
    # update values for next iteration
    prev_prev_num = prev_num 
    prev_num = current_num
}
```
<br>

## Task 3.4

Boxplot to show miles per gallon (in the variable mpg) as a function of the number of gears. Data taken from the `mtcars` dataset. 

```{r echo=T, message=F, warning=F}
library(tidyverse)
```


```{r mtcars}
ggplot(mtcars) + 
  geom_boxplot(aes(x = factor(gear), y = mpg, fill = factor(gear))) + 
  ylab("Miles per gallon (mpg)") + xlab("No. of gears") +
  theme_bw() + theme(legend.position = "none")
```
<br>
<br>

## Task 3.5

Performing a linear regression of car speed versus breaking distance using the `lm()` function. Data taken from the `cars` dataset.   

```{r}
# assign linear regression to a variable
reg_output = summary(lm(formula= dist~speed, data = cars))

# assign linear regression coefficients to variables
reg_coeffs = reg_output$coefficients

# extract and print slope values
slope = reg_coeffs["speed","Estimate"]
slope_se = reg_coeffs["speed","Std. Error"]

print(paste("The fitted slope of the line is", round(slope, digits = 2), "with a standard error of", round(slope_se, digits = 2)))

# extract and print intercept values
intercept = reg_coeffs["(Intercept)","Estimate"]
intercept_se = reg_coeffs["(Intercept)","Std. Error"]

print(paste("The intercept of the line is", round(intercept,digits = 2), "with a standard error of",  round(intercept_se,digits = 2)))
```
The units in the `cars` dataset are miles per hour (mph) for speed and feet (ft) for stopping distance (dist). Indentified in the `?cars` help documentation.
<br>
<br>

## Task 3.6

Plotting the data points from `cars` and the linear fit calculated in Task 3.5 with ggplot. 
```{r}
ggplot(cars, mapping = aes(x = speed, y = dist)) +
  geom_abline(intercept = intercept, slope = slope, colour = "black") + # these values are taken directly from task 3.5
  geom_point(colour = "dodgerblue4") +
  xlab("Speed (mph)") + ylab("Stopping distance (ft)") +
  theme_bw() 
```
<br>
<br>

## Task 3.7

Using linear regression to estimate the average reaction time for a driver to start braking (in seconds). Data taken from the `cars` dataset. To simplify matters it is assumed that once braking commences, braking distance is proportional to the square of the speed. 

As braking distance is proportional to the square of the speed, this can be modelled in a linear regression:
```{r}
cars_edit = cars %>% # Adding the square of the speed to the dataset
  mutate(speed_squared = speed^2)

reaction_reg = summary(lm(formula = dist ~ speed_squared, data = cars_edit))
```
<br>

The slope of this linear regression describes the relationship between stopping distance and speed squared. The y-intercept (where mph<sup>2</sup> is equal to 0) therefore describes the theoretical mean value of braking distance when speed is equal to 0, i.e. the portion of braking distance that is not dependent on speed: the reaction distance.    
```{r}
rxn_intercept = reaction_reg$coefficients["(Intercept)","Estimate"]
rxn_intercept_se = reaction_reg$coefficients["(Intercept)","Std. Error"]

print(paste("The intercept of the linear fit is", round(rxn_intercept,digits = 2), "ft with a standard error of",  round(rxn_intercept_se,digits = 2),"ft."))
```
<br>

As this value is taken from the linear fit of the data points in `cars`, we can assume that this average is the average distance a car going at mean average speed in the dataset travels before starting to brake (i.e. the reaction distance). 
```{r}
mean_speed = mean(cars_edit$speed)
print(paste("The mean average speed of the datapoints in cars is", round(mean_speed,digits = 2), "mph."))
```
<br>

We can now use the mean speed and reaction distance to estimate the average reaction time, by rearranging values into the equation "speed = distance/time". As our intercept is in feet, we will convert our mean speed from mph into feet per second, by multiplying our speed value by 1.47.
```{r}
rxn_time = rxn_intercept / (mean_speed*1.47) 
print(paste("The estimated average of reaction time in cars is", round(rxn_time,digits = 2), "seconds."))
```
<br>

However, due to the significant standard error of the intercept in the linear fit, the average distance before braking commences could be anywhere between the upper and lower limits of this standard error. This can be used to calculate the upper and lower limits of average reaction: 
```{r}
# Calculate upper limit
rxn_intercept_upper_lim = rxn_intercept + rxn_intercept_se
rxn_time_upper_lim = rxn_intercept_upper_lim / (mean_speed*1.47) 

print(paste("The upper limit of the estimated average of reaction time in cars is", round(rxn_time_upper_lim,digits = 2), "seconds."))

# Calculate lower limit
rxn_intercept_lower_lim = rxn_intercept - rxn_intercept_se
rxn_time_lower_lim = rxn_intercept_lower_lim / (mean_speed*1.47) 

print(paste("The lower limit of the estimated average of reaction time in cars is", round(rxn_time_lower_lim,digits = 2), "seconds."))

```
<br>

Mean average reaction time from 81 million online tests is 0.284 seconds (source: https://humanbenchmark.com/tests/reactiontime/statistics), which is within the range of our estimation. As our estimate was based on lots of assumptions and data from the 1920s, we can conclude that the use of linear regression here gives a reasonable estimation of reaction time. 
<br>

Finally, we can plot this linear regression on a graph of braking distance as a function of speed squared from `cars`: 

```{r}
rxn_slope = reaction_reg$coefficients["speed_squared","Estimate"]

ggplot(cars_edit, mapping = aes(x= speed_squared, y = dist)) +
  geom_point() +
  xlab("Speed (mph) squared") + ylab("Stopping distance (ft)") +
  geom_abline(intercept = rxn_intercept, slope = rxn_slope, colour = "red") +
  theme_bw() 
```

<!-- The slope of this linear regression describes the relationship between stopping distance and speed squared. As this value is the gradient of a straight line on a 2-dimensional graph, its units are in y/x, in this case ft/mph<sup>2</sup> (stopping distance/speed<sup>2</sup>). We can convert this value to seconds to estimate reaction time: -->
<!-- ```{r} -->
<!-- reaction_slope = reaction_reg$coefficients["speed_squared", "Estimate"] # slope in ft/mph^2 -->

<!-- # values for conversion -->
<!-- feet_per_mile = 5280 -->
<!-- seconds_per_hour = 60*60 -->

<!-- # convert value to seconds -->
<!-- # answer = (1 / (sqrt(reaction_slope / feet_per_mile)))/seconds_per_hour -->
<!-- answer = sqrt((reaction_slope / feet_per_mile) * seconds_per_hour) -->

<!-- print(paste("The average reaction time is", round(answer,digits = 3), "seconds.")) -->

<!-- ``` -->





