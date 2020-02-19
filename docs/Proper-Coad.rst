Proper Coadition
================

*This is the Co-Addition (or image stacking) tool, below describes the algorithm and how to use it*
The following was developed from `B.Zackay and E.O. Ofek. (2015a) <https://arxiv.org/abs/1512.06872>`_. and `B.Zackay and E.O. Ofek. (2015b) <https://arxiv.org/abs/1512.06879>`_.



The Maths
---------

To find proper image from co-addition, *R*:

.. math::
   
   \widehat{R} = \frac{\sum_j  \frac{F_j}{\sigma_j} \overline{{\widehat{P_j}}} \widehat{M_j}} {\sum_j \sqrt{\frac{F_j^2}{\sigma_j^2} |\widehat{P_j}|}}

The derrivation can be found As described in  `B.Zackay, et al. (2016) <http://iopscience.iop.org/article/10.3847/0004-637X/830/1/27/pdf>`_. 
