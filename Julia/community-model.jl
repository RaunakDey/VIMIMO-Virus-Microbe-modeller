using MAT
using Plots

# Load the .mat file
mat_data = matread("/Users/rdey33/Downloads/test/VIMIMO-Virus-Microbe-modeller/data/qpcr.mat")

# Access the data struct
data = mat_data["data"]

# Extract time data (assuming it's stored in data.xdata)
time = data["xdata"]

# Extract the 10 time series (assuming they are stored in data.ydata)
time_series = data["ydata"]

# Create subplots in a (5,2) grid
plot_arr = plot(layout = (5, 2), size = (800, 600))

# Plot each time series on a logarithmic scale
for i in 1:10
    plot!(plot_arr[i], time, time_series[:, i], seriestype = :scatter, title = "Series $i", xlabel = "Time", ylabel = "Value")
end

# Show the plots
display(plot_arr)


