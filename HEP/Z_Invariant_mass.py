# This script is using CMS open data to predict and calculate to invariant mass of the z bosons
# The mass can be calculated based in the transverse momentum of muons and both pseudorapidity and angles (direction of muon)

dataset = pd.read_csv('https://github.com/MohamedElashri/Physics/blob/master/HEP/data.csv')

# Browse fist few files of the database
dataset.head()

# Invariant mass formula
invariant_mass = np.sqrt(2*dataset.pt1*dataset.pt2*(np.cosh(dataset.eta1-dataset.eta2) - np.cos(dataset.phi1-dataset.phi2)))

# plot histogram of invariant mass
plt.hist(invariant_mass, bins=500)
plt.xlabel('Invariant mass [GeV]')
plt.ylabel('Number of events')
plt.title('The histogram of the invariant masses of the two muons \n') 
plt.show()

