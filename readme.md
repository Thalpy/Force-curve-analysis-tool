# Setup
## Run this in terminal in the directory of the project to create a virtual environment and install the required packages:

`python3 -m venv myenv`

`myenv\Scripts\activate`

`pip install -r requirements.txt`

# How to use:

Run GUI_v2.py

## Parameters

    binsize: This parameter defines the size of the bins used for averaging data points. A larger binsize will result in a smoother curve but may obscure finer details.

    dfit_win: This parameter is used in the remove_farfield_drift function. It specifies the window size for fitting a line to a portion of the data far from the contact point to remove any background drift in the deflection signal. It's essentially defining the range of data used for this linear fit.
    The approach phase window.

    dfit_off: Also used in the remove_farfield_drift function, this parameter defines the offset from the guessed contact point to the start of the window for drift fitting. It helps in avoiding the influence of contact-related changes while fitting the drift.

    cfit_min and cfit_max: These parameters define the minimum and maximum deflection values within which the script tries to fit a line to identify the contact point. This fitting is important for accurately determining where the tip makes contact with the surface.

    fitbin: This parameter is used in the binning process (bin_z_df function) and defines the number of data points per bin. It's similar to binsize but specifically for the procedure that aligns the contact region.

    cthresh: This is a threshold for deflection used as a 'first pass' guess to identify the contact point. It helps in determining the initial guess for where the tip makes contact with the surface.

    Moving the Orange Area Left and Right: The orange area likely represents a range of data points used for a specific calculation, possibly related to the background drift removal or finding the contact point. To move this area left or right, you would adjust the parameters that define the window of data points considered for this part of the analysis. This could be the dfit_win and dfit_off parameters in your script, which define the window and offset for fitting the background drift.
    dfit_win = size of window
    dfit_off = offset from guessed contact point (Higher means more to the left)


    Moving the Purple Area Up and Down: The purple area seems to represent data points after the contact point where the tip is interacting with the surface. Moving this area up or down would involve adjusting the baseline or the deflection offset. This could be done by changing the way the background drift is removed or by adjusting the deflection threshold (cthresh) for contact.
    cfit_min = minimum deflection value for fitting a line to identify the contact point
    cfit_max = maximum deflection value for fitting a line to identify the contact point

    Increasing the Dots in the Purple Area: The density of dots in the purple area corresponds to the number of data points in the region of the curve that has been plotted. To increase the number of dots, you would likely need to reduce the fitbin parameter, which controls the number of data points per bin when averaging to smooth out the curve. A smaller fitbin value would result in less averaging and more individual data points being plotted.
