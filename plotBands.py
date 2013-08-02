import matplotlib.pyplot as plt
import matplotlib.cm as cm
import parseProcar, surface

def PlotBands(procar):
    # Iterate over bands to collect data.
    for i in range(1, procar.Nb+1):
        # Iterate over k-points to build (kx, energy) pairs for this band.
        # Entries in surface are 0 (bulk states) or 1 (surface states).
        # surface controls the color of the markers.
        kx, energy, surface = [], [], []
        for k in procar.kPoints:
            kx.append(k.kx)
            b = k.Band(i)
            energy.append(b.energy)
            if hasattr(b, 'surface') and b.surface:
                print(str(k.kx) + " " + str(i))
                surface.append(0.0)
            else:
                surface.append(0.9)
        # plot band as thin green line (to clarify connection btwn markers)
        plt.plot(kx, energy, 'g-', linewidth=1)
        # plot band as markers: bulk -> green, surface -> black
        plt.scatter(kx, energy, c=surface, facecolors=surface, edgecolor='b',
                    marker='o', cmap=cm.gray, linewidth=1, s=50)

    # Display the plot
    plt.show()

if __name__ == "__main__":
    procar = None
    with open('PROCAR', 'r') as f:
        procar = parseProcar.PROCAR(f, nonCol=True, lmDecomposed=False, storeIds=True)
    surface.MarkSurfaceStates(procar, 1, 0.9999)
    PlotBands(procar)
