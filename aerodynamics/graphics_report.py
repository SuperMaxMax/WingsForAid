from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate

y_s_01 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
l_01 = [0.0,42.41358774253995,79.73903869827488,112.32130285337455,142.10028563953725,170.97220213107886,200.40901074718667,231.40459162532002,264.589459509799,300.32296992147974,338.78781700355415,380.03252464349146,424.0203508837592,470.63564894637267,519.7135457050697,571.0273412966508,624.31411441116,679.2488186389917,735.4733898148925,792.5560339585412,850.0293082552931,907.3267039662062,963.8363885519854,1018.7953248553503,1071.3699233586353,1120.4527546633997,1164.7909043805846,1202.4709436941243,1231.1049547935816,1245.0266204977872]

y_s_025 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
l_025 = [0.0,56.312622019433235,108.51322555840858,155.7344574223505,198.47818773566792,237.72830713839065,274.59101765385617,310.0831791312885,345.06381248654486,380.20497022383563,416.0133867812102,452.8418872894931,490.9255322013042,530.3892872672976,571.2806146451184,613.5644479642635,657.1527201817125,701.8857961487099,747.5633893899898,793.9106107940972,840.6154925523425,887.273050546035,933.4360947377601,978.5198872893131,1021.8769622651442,1062.612749475054,1099.703637805931,1131.5283305138917,1156.0372922319793,1168.1893963813288]

#y_s_03 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
#l_03 = [0.0,73.3739669331326,141.1485677821,201.8043650912512,255.60109761005063,303.50099353540884,346.71838900679205,386.4400962311308,423.71211290101667,459.38574000932704,494.1236433704433,528.4106537809151,562.5874061958019,596.8674060454746,631.3695179010533,666.1257387012402,701.108394086444,736.2264571542815,771.3502172460769,806.2957647015554,840.8512583411863,874.746622442297,907.6862682474436,939.2943983413496,969.1611545453572,996.7345906534391,1021.3920645858468,1042.1634483448734,1057.8327630392196,1065.4180004874206]
#
#y_s_04 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
#l_04 = [0.0,78.50732417267119,151.77596856094002,218.08800840071225,277.32511327357764,330.0841757441667,377.298858665795,419.98361438501996,459.10959829870126,495.5299602242258,529.9612109980359,562.9745580827163,595.0134112941727,626.401116444619,657.3660583824224,688.0472213390078,718.5185243242129,748.7865196400579,778.8136664036045,808.506028504065,837.7380191225341,866.3268749216525,894.062858746818,920.6619329695097,945.8075967140988,969.0562950326364,989.9014428415663,1007.5297106679571,1020.91144796385,1027.4558873300753]

y_s_05 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
l_05 = [0.0,67.28462083730803,131.1273027499622,190.2997175669709,244.60485047665185,294.3388150316905,340.0612344086236,382.4287349374665,422.11024278876295,459.72497840752084,495.8231801677597,530.866183851698,565.2351001008334,599.2237071043961,633.0575280456759,666.886047597925,700.8059045513264,734.846199543723,768.9952013529063,803.1744131879066,837.2711494104104,871.0949586305317,904.4208719749963,936.9139097380807,968.1916290521091,997.6768176164594,1024.695015629211,1048.0978690042175,1066.4009659446594,1075.7052229946105]

y_s_065 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
l_065 = [0.0,71.51244506368805,139.77680440422023,203.491927789483,262.28104881997405,316.23587815776835,365.72411738502603,411.24492322154094,453.35040597993685,492.5843377486056,529.4569288646406,564.4199160985862,597.8672919687629,630.1243133056694,661.4600894163608,692.0788126574524,722.1380465260318,751.7358993264454,780.9334865701676,809.7335680842533,838.1084702539811,865.9643656820256,893.1783278346965,919.5361828937558,944.7858057512449,968.5153448733305,990.2358832713523,1009.0688970488185,1023.8622273244541,1031.4547447974376]

y_s_075 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
l_075 = [0.0,73.83978797677747,144.52090590639995,210.71014612023512,271.9475686712481,328.2267362092237,379.81893905542074,427.13972393587227,470.6743713282601,510.9185772324053,548.3513824385533,583.4091501079965,616.4827216657685,647.9053061979054,677.9615286624384,706.8787140441698,734.8420044797621,761.9830847310004,788.3995477697229,814.1367849728624,839.2123545789914,863.5857980740374,887.1909792400297,909.8834002322303,931.4865530929301,951.6883601020007,970.1128200085103,986.0534665455056,998.5734594495229,1005.0191768189067]

#y_s_08 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
#l_08 = [0.0,90.59526369514606,176.62435616193955,256.11174074313334,328.3221854210822,393.17085681504165,450.97524780916586,502.27076087424615,547.7022409462317,587.944050273735,623.6558767913607,655.4510582405048,683.8846778017015,709.4440890270035,732.5519375888167,753.565133472679,772.7839184819734,790.4522477763566,806.7690474995493,821.8863706652107,835.9218190821113,848.9519411257869,861.0264773751412,872.1533670702045,882.3169133647881,891.4449897058397,899.4351518096545,906.0665837020398,911.0355677524425,913.45659095433]

#y_s_09 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
#l_09 = [0.0,92.61389954758415,180.74464632386253,262.3901361495301,336.7416854307493,403.6263292511418,463.273873692619,516.1429533840592,562.8167058948087,603.9252755841462,640.1002948268392,671.9425176909247,700.0065754514919,724.7904508474644,746.7348458323016,766.222001160293,783.5812798116679,799.0907661990067,812.9850483317899,825.4562446335502,836.6624857575508,846.7264094543708,855.7442888534572,863.7799050361848,870.8752945656272,877.03516654554,882.2413628954172,886.4087461025939,889.4041427939976,890.7916703672577]

y_s_1 = [4.763469717548333,4.756483678839986,4.735546053970155,4.700718256600087,4.652102442660222,4.589841210709242,4.514117183667927,4.425152473154669,4.323208027993823,4.208582868807848,4.081613210938289,3.9426714782682666,3.7921652108390496,3.6305358694648997,3.4582575408524394,3.275835547022658,3.083804963114358,2.8827290479165786,2.673197591733526,2.455825186427982,2.2312494227174846,2.000129020010897,1.763141894270883,1.520983169569587,1.2743631391699621,1.0240051821132474,0.7706436414236031,0.5150216701534963,0.25788905158775427,2.9167839712154564e-16]
l_1 = [0.0,78.56533257259713,154.11658377693246,225.26558155801393,291.41026824857045,352.3725211378336,408.2490212051496,459.2973290331567,505.86892231467675,548.3557708671702,587.1610918865947,622.6739716016017,655.2603894782277,685.2512827817793,712.944349318313,738.5973813327081,762.435799128222,784.6464245159945,805.3885336290278,824.7847801559046,842.9354602884671,859.903235660204,875.7318620116924,890.4188836750641,903.9419838075695,916.2043883008597,927.0749902319089,936.2462573583603,943.2895773318504,946.8529366581311]

u=0     #x-position of the center
v=0    #y-position of the center
a=max(y_s_01)     #radius on the x-axis
b=0.9162*max(l_05)    #radius on the y-axis
# print(np.pi*a*b/4)
# print(-integrate.simps(l_05, y_s_05))

t = np.linspace(0, 0.5*np.pi, 30)

plt.figure(figsize=(8, 6))
plt.plot(u+a*np.cos(t), v+b*np.sin(t), label = "Elliptical Lift Distribution", color = 'r', linestyle = 'dashed', linewidth=2)
plt.plot(y_s_01, l_01, label = r"$\lambda$ = 0.1", marker = '^', markevery=7, color = 'g', linewidth=1)
#plt.plot(y_s_025, l_025, label = r"$\lambda$ = 0.25", marker = 'v')
#plt.plot(y_s_03, l_03, label = r"$\lambda$ = 0.3")
#plt.plot(y_s_04, l_04, label = r"$\lambda$ = 0.4")
#plt.plot(y_s_05, l_05, label = r"$\lambda$ = 0.5", marker = '^')
#plt.plot(y_s_06, l_06, label = r"$\lambda$ = 0.6")
plt.plot(y_s_065, l_065, label = r"$\lambda$ = 0.65", marker = 's', markevery=7, color = 'k', linewidth=1)
#plt.plot(y_s_08, l_08, label = r"$\lambda$ = 0.8")
#plt.plot(y_s_09, l_09, label = r"$\lambda$ = 0.9")
plt.plot(y_s_1, l_1, label = r"$\lambda$ = 1.0",marker = 'v', markevery=7, color = 'b', linewidth=1)
plt.xlabel("Wing Span [m]")
plt.ylabel("Sectional Lift [N/m]")
plt.legend()
plt.savefig('spanwise_lift_changing_taper.pdf')
plt.show()

# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.8 indidence =  0.9852341903547764
# Span_eff =  0.9923262827208658 CL_wing =  0.38050604379776876 CL required for cruis =  0.380505722301332 CD_i =  0.005992624275161098
# 0.9923262827208659

# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.7 indidence =  0.9518558130081114
# Span_eff =  0.9938013748589707 CL_wing =  0.380503853793402 CL required for cruis =  0.380505722301332 CD_i =  0.005983660588053073
# 0.9938013748589711

# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.65 indidence =  0.9354402175917187
# Span_eff =  0.9939869778685633 CL_wing =  0.3805039868428525 CL required for cruis =  0.380505722301332 CD_i =  0.005982547468057825
# 0.9939869778685634

# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.6 indidence =  0.9193728317750671
# Span_eff =  0.99371666238576 CL_wing =  0.38050592790373716 CL required for cruis =  0.380505722301332 CD_i =  0.005984235922942348
# 0.9937166623857603

# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.5 indidence =  0.8888298754548696
# Span_eff =  0.9914544353872745 CL_wing =  0.3805057996157076 CL required for cruis =  0.380505722301332 CD_i =  0.005997886263046358
# 0.9914544353872745

# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.4 indidence =  0.8621172247318302
# Span_eff =  0.9860752145215685 CL_wing =  0.38050628495179845 CL required for cruis =  0.380505722301332 CD_i =  0.00603062121506182
# 0.9860752145215683

# =====================================================================
# current option is: AR =  6 taper ratio =  0.65 indidence =  1.201522096022886
# Span_eff =  0.9957311390225266 CL_wing =  0.3805047765660765 CL required for cruis =  0.380505722301332 CD_i =  0.007713953453594521
    
# =====================================================================
# current option is: AR =  7.75 taper ratio =  0.65 indidence =  0.9376538357615047
# Span_eff =  0.9939869688722845 CL_wing =  0.380505105078294 CL required for cruis =  0.380505722301332 CD_i =  0.0059825826856017335

# =====================================================================
# current option is: AR =  9.5 taper ratio =  0.65 indidence =  0.7703391230856338
# Span_eff =  0.9922470447625952 CL_wing =  0.3805038522500901 CL required for cruis =  0.380505722301332 CD_i =  0.0048890538842496795