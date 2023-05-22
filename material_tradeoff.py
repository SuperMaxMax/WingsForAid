import numpy as np 

### ALUMINUM ###
#1. 2024
#2. 2014
#3. 3003
#4. 5052
#5. 6061
#6. 7050
#7. 7068
#8. 7075

### STAINLESS STEEL ###
#9. 301
#10. 304
#11. 310
#12. 316
#13. 321
#14. 347
#15. 410
#16. 430

### STEEL ###
#17. 4130
#18. 4340

### NICKEL ###
#19. 200
#20. 400
#21. 600
#22. 625
#23. 718
#24. A286
#25. C276
#26. K500
#27. R405

### TITANIUM ###
#28. 6AL-4V
#29. Grade 2
#30. Grade 5

### COMPOSITES ###
#31. CFRP
#32. AFRP
#33. GFRP
#34. GLARE
#35. ARALL

### PLASTICS ###
#36. Polymethyl Methacrylate (PMMA)
#37. Thermosetting Polyimide
#38. Polyamide-imide (PAI)
#39. Polychlorotrifluoroethylene (PCTFE)
#40. Polytetrafluoroethylene (PTFE)
#41. Polyetheretherketon (PEEK)
#42. Polyphenylsulfon (PPSU)
#43. Ethyleen ChloorTriFluor Ethyleen (ECTFE)

### MICS ###
#44. Sitka spruce 
#45. Birch
#46. Ash
#47. Douglas fir
#48. Cupper
#49. Tungsten
#50. Magnesium

materials = np.array((),
                     (),
                     (),
                     ()
                     )

weights_fuel_tank = [0, 0.5 0.7, ]

score = 0
for i in len(materials):
    score_new = materials[i]*weights_fuel_tank
    if score_new > score: 
        score = score_new
        material_choice = materials[i]

print("yey you have selected", material_choice)