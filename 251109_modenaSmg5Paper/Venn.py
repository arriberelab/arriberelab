# Import necessary libraries
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Function to read a file and return a set of lines
def read_file(file_path):
    with open(file_path, 'r') as file:
        return set(file.read().splitlines())

# Read the files into sets (replace these with your actual file paths)
set1 = read_file('5m_targetList.txt')
set2 = read_file('5d_targetList.txt')
set3 = read_file('6d_targetList.txt')

# Create the Venn diagram
venn = venn3([set1, set2, set3], set_labels=('smg-5(D464A)', 'smg-5(Δ)', 'smg-6(Δ)'))

for i in ['100', '010', '001', '110', '101', '011', '111']:
    venn.get_patch_by_id(i).set_edgecolor('black')  # optional, adds black edges to patches
    venn.get_patch_by_id(i).set_linewidth(2)       # optional, to increase the line width of the edges

# For example, you can change labels or colors
venn.get_label_by_id('100').set_text('smg-5(D464A) only')
venn.get_label_by_id('010').set_text('smg-5(Δ) only')
venn.get_label_by_id('001').set_text('smg-6(Δ) only')

# Customize colors (optional)
venn.get_patch_by_id('100').set_color('red')
venn.get_patch_by_id('010').set_color('green')
venn.get_patch_by_id('001').set_color('blue')
#save plot
plt.savefig('venn_diagram_normal.svg', format='svg')
# Display the Venn diagram
plt.show()
