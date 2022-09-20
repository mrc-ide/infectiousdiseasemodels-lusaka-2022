## Example 3: Imperfect vaccines

Vaccination efficacy in the real world is not perfect. We are now going to model 'imperfect' vaccines that can either offer an 'all-or-nothing' protection or a 'leaky' protection.

An 'all-or-nothing' vaccine will be one where a proportion of those receiving the vaccine will receive the full protection effect, whereas the rest will have no protection at all.

On the other hand, a 'leaky' vaccine is one where all individuals receive some protection. However, they will 'leak' out of the \(Vaccinated\) compartment at a rate

$$\theta = \beta * \frac{I}{N} * (1 - ve)$$

where \(ve\) is the vaccine efficacy parameter, which is a proportion. 

You are provided with the code for an 'all-or-nothing' vaccine. Use this for part one of this example and then follow the instructions in the handout provided to help you code a 'leaky' vaccine in part two of the example.

