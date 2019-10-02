---
geometry: a4paper, margin=3cm
---

# Separating change in colony size into biological growth and partial mortality

Joshua S. Madin, Andrew H. Baird, Marissa Baskett, Sean R. Connolly, Maria Dornelas

#### Abstract

#### Introduction

- Growth is a fundamental process in all organisms. (will be elaborated)
- Growth in colonial organisms is complicated by the fact that biological growth (i.e., the addition of colony units) can be counteracted by the loss of units, known as partial mortality.
- When the loss of units is greater than the addition of units, colonies shrink.
- In most situations, these counteracting processes can be ignored. For instance, demographic modeling of colonial organisms is typically only concerned with changes in colony size through time, and less-so with how these changes come about.
- However, situations arise when it is important to separate the effects of biological growth and partial mortality because they are affected by different environmental and physiological variables. Therefore, understanding the impact of these variables on population dynamics requires that they be modelled separately.
- In this study, we developed a general model for separating change in colony size into biological growth and partial mortality components. We develop the model with 11 species of reef building corals that exhibit five distinctly different growth morphologies. Such a model will help understand and quantify how environmental changes will influence overall colony growth and subsequently population persistence.

#### Methods

We estimated yearly changes in planar area for 11 species of scleractinian corals, across five morphological groups: tabular (Acropora cytherea and A. hyacinthus), arborescent (A. intermedia and A. robusta), corymbose (A. spathulata, A. millepora and A. nasuta), digitate (A. cf digitifera and A. humilis) and massive (Goniastrea pectinata and G. retiformis). We tagged 30 colonies of each species distributed along a 500 m by 10 m band of reef crest at Trimodal reef, between South and Palfrey islands (14.6998398 S, 145.4486748 E), Lizard Island, Australia. Each tagged colony was photographed from above with a two-dimensional scale plate placed level with the surface of the colony. The angle of the camera was horizontal, and the distance from the colony was such that the entire colony was visible in the photograph. Colonies were photographed yearly in November 2009, 2010, 2011, 2012 and 2013. The images were corrected for barrel distortion, and the scale and outline of each colony were digitized in IMAGEJ for estimation of planar area. Every year, dead or missing colonies were replaced in order to maintain approximately 30 colonies per species. To minimize the effect of observation error, colonies were photographed twice to three times independently every year. All photographs were digitized twice independently and estimates of area were averaged.

Biological growth occurs at the periphery of colonies, and so we calculated the changes in the radius of colonies each year. Radii were calculated from planar areas by assuming colonies were circular. We also calculated the circularity of colonies, which is the ratio of colony outline lengths and the perimeter of a circle with the same planar area as the colony. These estimates of circularity were used to test for biases in our analysis below.

We assumed that the upper bound of added radii as a function of colony planar area (i.e., composed from the fastest growing individuals in the population) is reflective of typical biological growth in the absence of significant partial mortality. We use quantile regression to detect this upper bound. We use the *rq* function in the package *quantreg* (ref) in R (ref) to calculated the upper 97.5% quantile and determined if the slope of the quantile is statistically indistinguishable from zero (i.e., added radius is size independent). Partial mortality is bounded between biological growth (i.e., no partial mortality) and whole-colony death (i.e., partial mortality that is equal to colony size). We used the beta function to model partial mortality as the proportion of planar area lost relative to the mean potential biological growth per year.

#### Results

The relationship between radial growth and colony size (planar area) yielded relationships with constant maximum radial growth for ten out of the 11 species (Table 1; Fig. 1). The exception, *Acropora spathulata*, tended to add more radius per year when larger. Circularity ranged broadly for some species (Fig. 2); however, we found no strong associations between circularity and model residuals, suggesting that our results are robust to calculating added radius based on colonies being circular. Biological growth was greatest in the tabular species (Acropora hyacnthus and A. cytherea) that tended to add approximatley 10 cm in radius per year. Added radius was lowest for the massive species (Goniarstrea retiformis and G. pectinata) at about 2 cm per year.

Fits of the beta distirbution to the the proportions of partial mortality experience by colonies tended to be modal for species; however, the peakedness and skewless of this peak differed among taxa. For example, partial mortality tended to 


#### Discussion

-


Caveats:
- Using quantile regression to estimate biological growth will not capture true biological growth, because we do not know the true amount of biological growth and partial mortality for each individual. However, in the absence of experiments that minimise partial mortality and carefully monitor growth through time, we believe this is the best we can do.
- Coral colonies are not circular. Assuming they are will underestimate added radius in less circular colonies.


Table 1: Regression of upper 0.975 quantiles, illustrating that the slopes were on the whole indistinguishable for zero. *Acropora spathulata* had a significant slope.

Species |Value | Std. Error | t value | Pr(>)
--- | --- | --- | --- | ---
*Acropora hyacinthus* | -0.006 | 0.015 | -0.425 | 0.672
*Acropora cytherea* | 0.059 | 0.032 | 1.840 | 0.070 .
*Acropora intermedia* | 0.019 | 0.011 | 1.666 | 0.102
*Acropora robusta* | 0.028 | 0.018 | 1.507 | 0.135
*Acropora nasuta* | -0.002 | 0.009 | -0.177 | 0.860
*Acropora spathulata* | 0.030 | 0.011 | 2.656 | 0.009 **
*Acropora millepora* | 0.000 | 0.011 | 0.033 | 0.974
*Acropora fat dig* | -0.007 | 0.013 | -0.512 | 0.610
*Acropora humilis* | 0.010 | 0.012 | 0.861 | 0.392
*Goniastrea retiformis* | 0.004 | 0.004 | 0.886 | 0.378
*Goniastrea pectinata* | -0.018 | 0.011 | -1.563 | 0.121

Decrease in radius (partial mortality) was looked at proportionally, which assumes the same patterns of partial mortality are size independent because they are bounded.

![The estimated proportion of colony partial mortality following expected yearly growth and beta distribution model fit (solid curves).](output/beta_fit.png)

The different growth forms have different shaped partial mortality distributions, which are essential colony area multipliers (Figure 2); i.e., they range from almost full colony mortality at the zero end, through to no partial mortality at the one end. For example, *Acropora hyacinthus* colonies tend to peak around 0.5, meaning that on average a colony loses half its area to partial mortality. However, note that this estimate applies *after* yearly growth, and so this doesn't necessarily mean overall growth is negative, which is especially the case for smaller colonies. Other growth forms, like massives, tend to lose far less colony area on average. I fit a beta distribution to these partial mortality plots, because they were bounded by zero and one.

From the constant radial growth increment (modeled as normal distribution with the bootstrapped standard deviation) and the partial mortality multiplier (model the a beta distribution), the growth function is below:

    g_func <- function(x, y, qr0, shape1, shape2) {
      r0 <- r_func(x)

      # Actual radial growth
      r1 <- r0 + qr0

      # Probability of y?
      return(dbeta(y / a_func(r1), shape1, shape2))
    }

This give the probability of a colony of size *x* growth to size *y* in a year. I'm not sure how to express this as an equation.

Below are the resulting probability heat map with the data superimposed (Figure 3).

![Planar area of colonies at t+1 as a function of the area at time t. Dashed line in the zero growth line. Shading represents estimated probability of early growth (t+1) as a function of colony size at t (green is highest probability).](output/growth_relative.png)


![Increase in mean colony radius as a function of colony planar area. Upper dashed line is the 0.975 quantile regression with no slope parameter, because generally the slope was not discernible from 0; see Table. Shaded area is the 95% confidence interval for the quantile. The bottom dashed line represents the maximum amount of radius that can be lost before whole-colony mortality.](output/added_radius.png)
