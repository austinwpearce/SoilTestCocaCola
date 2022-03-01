# This testing file includes a soil test correlation dataset from the agridat
# package as an example for testing the lin_plateau, quad_plateau, mitscherlich,
# and ALCC functions.

cotton_k <- agridat::cate.potassium
corn_n <- agridat::engelstad.nitro # QP
yield_monitor <- agridat::gartner.corn
yield_monitor <- agridat::lasrosas.corn
agridat::sinclair.clover


cotton <- tibble(stk = agridat::cate.potassium$potassium,
                 ry = agridat::cate.potassium$yield)