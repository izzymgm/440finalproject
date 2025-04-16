#!/usr/bin/env python3
"""
plot_nod_vs_oxygen.py

Creates a scatter plot of the number of nod genes detected versus oxygen concentration,
fits a linear regression line, computes the R² score, and displays an anoxic cutoff.
"""

# Import necessary libraries
import os                                   
import matplotlib.pyplot as plt             
import numpy as np                         
import pandas as pd                        
from sklearn.metrics import r2_score        


def main():
    # -------------------------------------------------------------------------
    # 1. Load your data using a relative path (place df_all_o2readings.csv in same folder as script
    # -------------------------------------------------------------------------

    df_ALL = pd.read_csv("df_all_o2readings.csv")

    # -------------------------------------------------------------------------
    # 2. Extract the variables for plotting
    # -------------------------------------------------------------------------
    x = df_ALL["# of nod genes"]   # number of nod genes detected in each sample
    y = df_ALL["Oxygen"]           # corresponding oxygen concentration measurements

    # -------------------------------------------------------------------------
    # 3. Fit a simple linear regression
    # -------------------------------------------------------------------------
    slope, intercept = np.polyfit(x, y, 1)  # returns [slope, intercept]
    y_pred = slope * x + intercept         # predicted y values using the fitted line
    r_squared = r2_score(y, y_pred)        # coefficient of determination (R²)

    # -------------------------------------------------------------------------
    # 4. Configure global plotting parameters
    # -------------------------------------------------------------------------
    plt.rcParams.update({'font.size': 16})  # increase base font size for all text

    # -------------------------------------------------------------------------
    # 5. Create the figure and scatter plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 7))
    plt.scatter(
        x,
        y,
        color="blue",            # marker color
        alpha=0.4,                 # marker transparency
        s=80,                      # marker size
        label="environmental samples"
    )

    # -------------------------------------------------------------------------
    # 6. Plot the fitted regression line
    # -------------------------------------------------------------------------
    plt.plot(
        x,
        y_pred,
        color="green",           # line color
        linestyle="--",          # dashed line
        label=f"Linear fit ($R^2$ = {r_squared:.2f})"
    )

    # -------------------------------------------------------------------------
    # 7. Add an anoxic cutoff line at Oxygen = 20
    # -------------------------------------------------------------------------
    plt.axhline(
        y=20,
        color="red",
        linestyle="-.",
        label="anoxic cutoff"
    )

    # -------------------------------------------------------------------------
    # 8. Label axes and title
    # -------------------------------------------------------------------------
    plt.xlabel("# of nod genes detected in sample")
    plt.ylabel("Oxygen concentration")
    plt.title("Scatter Plot of # of nod genes vs Oxygen")

    # -------------------------------------------------------------------------
    # 9. Add legend and display
    # -------------------------------------------------------------------------
    plt.legend()
    plt.tight_layout()  # adjust spacing to prevent clipping of labels
    plt.show()


if __name__ == "__main__":
    main()
