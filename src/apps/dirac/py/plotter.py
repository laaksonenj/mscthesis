import argparse
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import Slider, CheckButtons, TextBox

parser = argparse.ArgumentParser()
parser.add_argument('--plot-data', nargs='+', required=True)
parser.add_argument('--labels', nargs='+')
parser.add_argument('--mesh-file')
parser.add_argument('--dirac-point', nargs='+')
args = parser.parse_args()

plot_data = []
for filename in args.plot_data:
    with open(filename) as f:
        plot_data.append(np.loadtxt(f))
labels = args.labels if args.labels is not None else []

assert len(labels) == 0 or len(labels) == len(plot_data)
assert not args.dirac_point or len(args.dirac_point) == 2

max_p = max([np.shape(arr)[0] for arr in plot_data])

# C1 * p^(-1) * (sqrt(log(p + 1)) + 1)
def thesis_bound(t, C1):
    return C1 * (1 / t) * (np.sqrt(np.log(t + 1)) + 1)

# C2 * p^(-k)
def algebraic_bound(t, C2, k):
    return C2 * np.power(t, -k)

fig = plt.figure(figsize=(9, 4.8))
fig.subplots_adjust(left=0.4)
ax = fig.subplots()

init_C1 = 0.25
init_C2 = 0.05
init_k = 1

t = np.linspace(0.9, max_p + 1, 100)
line_thesis_bound, = ax.plot(t, thesis_bound(t, init_C1), label=r'$C_1 p^{-1} ( \sqrt{\log (p + 1)} + 1 )$')
line_algebraic_bound, = ax.plot(t, algebraic_bound(t, init_C2, init_k), label=r'$C_2 p^{-k}$')

# Controls for thesis_bound
thesis_bound_checkbutton_ax = fig.add_axes([0.05, 0.85, 0.2, 0.05])
thesis_bound_checkbutton = CheckButtons(
    ax=thesis_bound_checkbutton_ax,
    labels=['Show thesis bound'],
    actives=[True]
)

C1_slider_ax = fig.add_axes([0.05, 0.8, 0.2, 0.03])
C1_slider = Slider(
    ax=C1_slider_ax,
    label='$C_1$',
    valmin=0.01,
    valmax=1,
    valinit=init_C1
)
C1_slider.valtext.set_visible(False)
C1_slider.vline.set_linewidth(0)

C1_textbox_ax = fig.add_axes([0.255, 0.8, 0.05, 0.03])
C1_textbox = TextBox(
    ax=C1_textbox_ax,
    label='',
    initial=str(init_C1),
    textalignment='center'
)

def thesis_bound_checkbutton_on_clicked(label):
    if line_thesis_bound.get_visible():
        line_thesis_bound.set_visible(False)
        line_thesis_bound.set_label('_' + line_thesis_bound.get_label())
    else:
        line_thesis_bound.set_visible(True)
        line_thesis_bound.set_label(line_thesis_bound.get_label()[1:])
    ax.legend()
    fig.canvas.draw_idle()

def C1_slider_on_changed(val):
    line_thesis_bound.set_ydata(thesis_bound(t, val))
    display_val = f'{val:.3f}'
    C1_textbox.set_val(display_val)
    fig.canvas.draw_idle()

def C1_textbox_on_submit(val):
    C1_slider.set_val(float(val))

thesis_bound_checkbutton.on_clicked(thesis_bound_checkbutton_on_clicked)
C1_slider.on_changed(C1_slider_on_changed)
C1_textbox.on_submit(C1_textbox_on_submit)

# Controls for algebraic bound
algebraic_bound_checkbutton_ax = fig.add_axes([0.05, 0.7, 0.25, 0.05])
algebraic_bound_checkbutton = CheckButtons(
    ax=algebraic_bound_checkbutton_ax,
    labels=['Show algebraic bound'],
    actives=[True]
)

C2_slider_ax = fig.add_axes([0.05, 0.65, 0.2, 0.03])
C2_slider = Slider(
    ax=C2_slider_ax,
    label='$C_2$',
    valmin=0.01,
    valmax=1,
    valinit=init_C2
)
C2_slider.valtext.set_visible(False)
C2_slider.vline.set_linewidth(0)

C2_textbox_ax = fig.add_axes([0.255, 0.65, 0.05, 0.03])
C2_textbox = TextBox(
    ax=C2_textbox_ax,
    label='',
    initial=str(init_C2),
    textalignment='center'
)

k_slider_ax = fig.add_axes([0.05, 0.6, 0.2, 0.03])
k_slider = Slider(
    ax=k_slider_ax,
    label='$k$',
    valmin=0.1,
    valmax=10,
    valinit=init_k
)
k_slider.valtext.set_visible(False)
k_slider.vline.set_linewidth(0)

k_textbox_ax = fig.add_axes([0.255, 0.6, 0.05, 0.03])
k_textbox = TextBox(
    ax=k_textbox_ax,
    label='',
    initial=str(init_k),
    textalignment='center'
)

def algebraic_bound_checkbutton_on_clicked(label):
    if line_algebraic_bound.get_visible():
        line_algebraic_bound.set_visible(False)
        line_algebraic_bound.set_label('_' + line_algebraic_bound.get_label())
    else:
        line_algebraic_bound.set_visible(True)
        line_algebraic_bound.set_label(line_algebraic_bound.get_label()[1:])
    ax.legend()
    fig.canvas.draw_idle()

def C2_slider_on_changed(val):
    line_algebraic_bound.set_ydata(algebraic_bound(t, val, k_slider.val))
    display_val = f'{val:.3f}'
    C2_textbox.set_val(display_val)
    fig.canvas.draw_idle()

def C2_textbox_on_submit(val):
    C2_slider.set_val(float(val))

def k_slider_on_changed(val):
    line_algebraic_bound.set_ydata(algebraic_bound(t, C2_slider.val, val))
    display_val = f'{val:.3f}'
    k_textbox.set_val(display_val)
    fig.canvas.draw_idle()

def k_textbox_on_submit(val):
    k_slider.set_val(float(val))

algebraic_bound_checkbutton.on_clicked(algebraic_bound_checkbutton_on_clicked)
C2_slider.on_changed(C2_slider_on_changed)
C2_textbox.on_submit(C2_textbox_on_submit)
k_slider.on_changed(k_slider_on_changed)
k_textbox.on_submit(k_textbox_on_submit)

# Axis scaling
x_axis_log_scale_checkbutton_ax = fig.add_axes([0.05, 0.5, 0.12, 0.05])
x_axis_log_scale_checkbutton = CheckButtons(
    ax=x_axis_log_scale_checkbutton_ax,
    labels=['log x-axis']
)

def x_axis_log_scale_checkbutton_on_clicked(label):
    if x_axis_log_scale_checkbutton.get_status()[0]:
        ax.set_xscale('log')
    else:
        ax.set_xscale('linear')
    fig.canvas.draw_idle()

x_axis_log_scale_checkbutton.on_clicked(x_axis_log_scale_checkbutton_on_clicked)

y_axis_log_scale_checkbutton_ax = fig.add_axes([0.19, 0.5, 0.12, 0.05])
y_axis_log_scale_checkbutton = CheckButtons(
    ax=y_axis_log_scale_checkbutton_ax,
    labels=['log y-axis']
)

def y_axis_log_scale_checkbutton_on_clicked(label):
    if y_axis_log_scale_checkbutton.get_status()[0]:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    fig.canvas.draw_idle()

y_axis_log_scale_checkbutton.on_clicked(y_axis_log_scale_checkbutton_on_clicked)

# Plot data
for data, label in zip(plot_data, labels):
    ax.plot(data[:,0], data[:,1], label=label, linestyle='-', marker='s', markerfacecolor='none')

# Mesh
def rational(val):
    try:
        return float(val)
    except ValueError:
        num, den = val.split('/')
        return float(num) / float(den)

if args.mesh_file:
    mesh_ax = fig.add_axes([0.05, 0.11, 0.3, 0.3])
    mesh_ax.set_aspect('equal', 'box')
    mesh_ax.set_title('Mesh')
    mesh_ax.set_xlim([-1.01, 1.01])
    mesh_ax.set_ylim([-1.01, 1.01])
    with open(args.mesh_file) as f:
        lines = [line.rstrip() for line in f]
    nodes = []
    elements = []
    for line in lines:
        if line[0] == 'n':
            x, y = line[2:].split()
            nodes.append(np.array([rational(x), rational(y)]))
        if line[0] == 'e':
            indices = line[2:].split()
            elements.append(np.array([int(x) for x in indices]))
    for element in elements:
        num = len(element)
        for i in range(num):
            n1 = nodes[element[i]]
            n2 = nodes[element[(i + 1) % num]]
            x = np.array([n1[0], n2[0]])
            y = np.array([n1[1], n2[1]])
            mesh_ax.plot(x, y, color='black')
    if args.dirac_point:
        dirac_point_x = rational(args.dirac_point[0])
        dirac_point_y = rational(args.dirac_point[1])
        mesh_ax.plot(dirac_point_x, dirac_point_y, 'ro')

ax.legend()
plt.show()
