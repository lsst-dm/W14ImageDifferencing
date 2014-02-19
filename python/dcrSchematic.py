import matplotlib.pyplot as plt
import numpy as np

def deg2rad(deg): return deg * np.pi / 180

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
ax.set_theta_zero_location('N')
plt.setp(ax.get_yticklabels(), visible=False) # turn off radial labels
plt.plot((0,0), (0,1), "k-", linewidth=3)
plt.plot((0,deg2rad(270)), (0,1), "k-", linewidth=3)
ax.set_xticklabels(("Y", "$45^\circ$", "$90^\circ$", "$135^\circ$", "$180^\circ$", "$225^\circ$", "X", "$315^\circ$"), weight="bold", fontsize=20)

# Pole
coord = (deg2rad(251.06), 0.9)
ax.annotate("", xytext=(0.0,0.0), xy=coord,
            arrowprops=dict(facecolor="red", width=4, headwidth=10))
ax.text(coord[0], coord[1], "P", color="red", fontsize=18, weight="bold", horizontalalignment="left", verticalalignment="center")

# DCR
coord = (deg2rad(145.71), 0.8)
ax.annotate("", xytext=(0.0,0.0), xy=coord,
            arrowprops=dict(facecolor="blue", width=3, headwidth=8))
ax.text(coord[0], coord[1], "A", color="blue", fontsize=18, weight="bold", horizontalalignment="center", verticalalignment="top")

coord = (deg2rad(152.27), 0.8)
ax.annotate("", xytext=(0.0,0.0), xy=coord,
            arrowprops=dict(facecolor="blue", width=3, headwidth=8))
ax.text(coord[0], coord[1], "B", color="blue", fontsize=18, weight="bold", horizontalalignment="center", verticalalignment="top")

coord = (deg2rad(349.86), 0.8)
ax.annotate("", xytext=(0.0,0.0), xy=coord,
            arrowprops=dict(facecolor="blue", width=3, headwidth=8))
ax.text(coord[0], coord[1], "D", color="blue", fontsize=18, weight="bold", horizontalalignment="center", verticalalignment="bottom")

coord = (deg2rad(356.41), 0.8)
ax.annotate("", xytext=(0.0,0.0), xy=coord,
            arrowprops=dict(facecolor="blue", width=3, headwidth=8))
ax.text(coord[0], coord[1], "E", color="blue", fontsize=18, weight="bold", horizontalalignment="center", verticalalignment="bottom")

#import matplotlib.patheffects as PathEffects
#t.set_path_effects([PathEffects.Stroke(linewidth=2, foreground="k"), PathEffects.Normal()])

# DCR angles
plt.plot((deg2rad(np.linspace(145.71, 251.06, 100))), 0.6 * np.ones(100), "g-",  linewidth=3)
plt.plot((deg2rad(np.linspace(251.06, 356.41, 100))), 0.6 * np.ones(100), "g--", linewidth=3)

plt.plot((deg2rad(np.linspace(152.27, 251.06, 100))), 0.4 * np.ones(100), "m-",  linewidth=3)
plt.plot((deg2rad(np.linspace(251.06, 349.86, 100))), 0.4 * np.ones(100), "m--", linewidth=3)

plt.show()
