## Debugging code 06/11/2024

Here I am trying to isolate the issue!

# Report 1: 06/11/2024. git commit report 1

There is a problem with the diffusion rate, the higher the diffusion then the hybrid method seems to exhibit higher profile which looks almost periodic. There is no problem with no production, the profiles are very similar. There is not PDE in this case. 


# Report 2: git commit report 2

Have set a low particle number per compartment, a conversion threshold. The solution still looks okay, seems to be behaving exactly as one would expect. May now increase the compartment numbers to see if thats where a problem occurs


# Report 3: git commit report 3

Raised the number of compartments to 12. The combined solution is only slightly higher than the analytic solution. Problem is beginning to emerge, must find the instance in which this breaks.

# Report 4: git commit report 4

I have decreased gamma to see if this a problem, I have also increased the number of repeats. This looks far better. No issues, next I am going to increase diffusion and reduce compartments 

# Report 5: git commit report 5

First decreased compartments from 12-->8. And increased diffusion from 10e-2--> 10e-1. Number of repeats from 20--> 10. Still no issues. Try to increase gamma. and decrease diffusion once again. 


# Report 6: git commit report 6

 