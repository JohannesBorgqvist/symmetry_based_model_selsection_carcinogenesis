import symmetry_toolbox # Home-made



# Calculate the midpoints
midpoint_myeloma = ((85-25)/(2))
midpoint_colon = ((85-12)/(2))
midpoint_CML = ((85-10)/(2))
# Calculate the factors
n =2
age = 85
# =================================================================================
# =================================================================================
# CALCULATE THE TRANSFORMATION SCALES
# =================================================================================
# What epsilon should we use to transform the age 85 years to 170 years? Here, comes
# the answer in the case of our two models.

# The PLM
epsilon_scale_PLM = symmetry_toolbox.PLM_transformation_scale(n)
# The IM-III
# Myeloma data
epsilon_scale_IM_III_myeloma = symmetry_toolbox.IM_III_transformation_scale(age,n,0.13,13.71)
# Colon data
epsilon_scale_IM_III_colon = symmetry_toolbox.IM_III_transformation_scale(age,n,0.097,33.24)
# CML data
epsilon_scale_IM_III_CML = symmetry_toolbox.IM_III_transformation_scale(age,n,0.064,-12.85)

print("\n\tEPSILON SCALES")
print("\n\t\tPLM\t\tn\t=\t%d:\t%0.12f"%(n,epsilon_scale_PLM))
print("\t\tIM-III myeloma\t(n,age)\t=\t(%0.3f,%d):\t%0.12f"%(n,age,epsilon_scale_IM_III_myeloma))
print("\t\tIM-III colon\t(n,age)\t=\t(%0.3f,%d):\t%0.12f"%(n,age,epsilon_scale_IM_III_colon))
print("\t\tIM-III CML\t(n,age)\t=\t(%0.3f,%d):\t%0.12f\n\n"%(n,age,epsilon_scale_IM_III_CML))
print("\n\tMIDPOINTS")
print("\t\tMidpoint myeloma\t=\t%0.3f"%(midpoint_myeloma))
print("\t\tMidpoint colon\t\t=\t%0.3f"%(midpoint_colon))
print("\t\tMidpoint CML\t\t=\t%0.3f"%(midpoint_CML))
