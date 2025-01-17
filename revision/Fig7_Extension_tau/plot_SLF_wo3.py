import numpy as np
import matplotlib.pyplot as plt

# Data from the table
a0 = np.array([-1, -0.5, 0, 0.5, 1])

folder_path = "revision/Fig7 Extension/wo3_SBD_NL_ke"
file_name = "ke_simu_SLF_choice.txt"
file_path = f"{folder_path}/{file_name}"
SBD_ke=np.loadtxt(file_path)

folder_path = "revision/Fig7 Extension/wo3_SBD_NL"
file_name = "simu_SLF_choice.txt"
file_path = f"{folder_path}/{file_name}"
SBD=np.loadtxt(file_path)

folder_path = "revision/Fig7 Extension/wo3_DLP_NL"
file_name = "simu_SLF_choice.txt"
file_path = f"{folder_path}/{file_name}"
DLP=np.loadtxt(file_path)


# Plot lines with dots
plt.figure(figsize=(8, 6))
plt.plot(a0, SBD_ke[0], '-o', label='SBD_ke')
plt.plot(a0, SBD[0], '-^', label='SBD')
plt.plot(a0, DLP[0], '-s', label='DLP')

# Add labels and title
plt.xlabel('Demand')
plt.ylabel('Revenue')
# plt.title('Increase of Relative Expected Revenue by selling Extra Seats over Demand Ratio (HOMOG)')
plt.xticks(a0, ['Strong, a0=-1', '', '', '', 'Weak, a0=1'])

# Add legend
plt.legend(loc='best')

# Adjust axes limits if necessary
plt.xlim([min(a0)-0.3, max(a0)+0.3])
plt.ylim([min(np.concatenate([SBD_ke[0], SBD[0], DLP[0]]))-1, 
          max(np.concatenate([SBD_ke[0], SBD[0], DLP[0]]))+1])

# Save the figure
#plt.savefig('Fig5.eps', format='eps')



# Plot lines with dots
plt.figure(figsize=(8, 6))
plt.plot(a0, SBD_ke[1], '-o', label='SBD_ke')
plt.plot(a0, SBD[1], '-^', label='SBD')
plt.plot(a0, DLP[1], '-s', label='DLP')

# Add labels and title
plt.xlabel('Demand')
plt.ylabel('SLF')
# plt.title('Increase of Relative Expected Revenue by selling Extra Seats over Demand Ratio (HOMOG)')
plt.xticks(a0, ['Strong, a0=-1', '', '', '', 'Weak, a0=1'])

# Add legend
plt.legend(loc='best')

# Adjust axes limits if necessary
plt.xlim([min(a0)-0.3, max(a0)+0.3])
plt.ylim([min(np.concatenate([SBD_ke[1], SBD[1], DLP[1]]))-1, 
          max(np.concatenate([SBD_ke[1], SBD[1], DLP[1]]))+1])

# Save the figure
#plt.savefig('Fig5.eps', format='eps')


# Plot lines with dots
plt.figure(figsize=(8, 6))
plt.plot(a0, SBD_ke[2], '-o', label='SBD_ke')
plt.plot(a0, SBD[2], '-^', label='SBD')
plt.plot(a0, DLP[2], '-s', label='DLP')

# Add labels and title
plt.xlabel('Demand')
plt.ylabel('SLF')
# plt.title('Increase of Relative Expected Revenue by selling Extra Seats over Demand Ratio (HOMOG)')
plt.xticks(a0, ['Strong, a0=-1', '', '', '', 'Weak, a0=1'])

# Add legend
plt.legend(loc='best')

# Adjust axes limits if necessary
plt.xlim([min(a0)-0.3, max(a0)+0.3])
plt.ylim([min(np.concatenate([SBD_ke[2], SBD[2], DLP[2]]))-1, 
          max(np.concatenate([SBD_ke[2], SBD[2], DLP[2]]))+1])

# Save the figure
#plt.savefig('Fig5.eps', format='eps')
plt.show()
