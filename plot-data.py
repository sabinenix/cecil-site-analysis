import geopandas as gpd
import pandas as pd
import yaml
import json
import os
import matplotlib.pyplot as plt
import seaborn as sns
from shapely import from_wkt
import cecil

# Import helper functions
from utils import show_linked_reprojections
from utils import get_reprojection_details_by_id
from utils import get_data_request_details_by_id
from utils import query_all_data_analytics
from utils import query_all_data_raw
from utils import df_to_gdf



def plot_all_vars(gdf, var_dict, year, provider_name, output_folder):
    """
    Plot all variables for one year given a dictionary of variable names and units.
    
    Args:
    gdf: GeoDataFrame created from the DataFrame of a single dataset pulled from Snowflake with 
    a shapely geometry column called 'geometry' created from the original column 'pixel_boundary'.
    var_dict: A dictionary with all the variables (keys) and units (values) for the dataset in gdf.
    year: The year for which to plot all the variables.
    provider_name: The name of the provider of the dataset. Note this is only used for labeling and naming
    the plots so they are more easily understood by human users.
    output_folder: The folder where the output plots should be saved.
    
    Outputs: A series of plots saved to output_folder.
    """
    # Check output folder exists and create it if it does not.
    os.makedirs(output_folder, exist_ok=True)

    # Subset data to just the year of interest.
    if 'year' in gdf.columns:
        gdf_year = gdf[gdf['year'] == year]
    elif 'date' in gdf.columns:
        gdf['date'] = pd.to_datetime(gdf['date'])
        gdf_year = gdf[gdf['date'].dt.year == year]
        gdf_year = gdf_year[gdf_year['date'] == gdf_year['date'].min()]
    else:
        raise ValueError("Neither 'year' nor 'date' column exists in the GeoDataFrame.")

    # Loop through variables and plot each
    for var in var_dict.keys():
        fig, ax = plt.subplots(figsize=(8,8))
        plot = gdf_year.plot(
            column = var, 
            ax=ax,
            markersize = 0.5,
            legend = True,
            cmap = "viridis"
            )
        
        # Pull colorbar units from var_dict.
        colorbar = plot.get_figure().axes[-1]
        colorbar.set_ylabel(f"{var_dict[var]}")

        # Set the plot options.
        ax.set_title(f"{provider_name} - {var}")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.tight_layout()

        # Save the plot to output folder.
        output_file = os.path.join(output_folder, f"{provider_name}_{var}.png")
        plt.savefig(output_file, dpi=300)
        print(f'Saved figure to: {output_file}')
        plt.close()


def plot_all_years(gdf, var_dict, var_of_interest, provider_name, output_folder):
    """
    Plot all years of a single variable.
    
    Args:
    gdf: GeoDataFrame created from the DataFrame of a single dataset pulled from Snowflake with 
    a shapely geometry column called 'geometry' created from the original column 'pixel_boundary'.
    var_dict: A dictionary with all the variables (keys) and units (values) for the dataset in gdf.
    var_of_interest: The variable for which to plot all the years.
    provider_name: The name of the provider of the dataset. Note this is only used for labeling and naming
    the plots so they are more easily understood by human users.
    output_folder: The folder where the output plots should be saved.
    
    Outputs: A series of plots saved to output_folder.
    """

    # Check output folder exists and create it if it does not.
    os.makedirs(output_folder, exist_ok=True)

    # Get the min and max values across timeseries to fix color ramp.
    max_val = gdf[var_of_interest].max()
    min_val = gdf[var_of_interest].min()
    print(f"Global range for {var_of_interest}: {min_val} to {max_val}")

    # Determine if 'year' or 'time' column exists.
    if 'year' in gdf.columns:
        group_col = 'year'
    elif 'date' in gdf.columns:
        gdf['date'] = pd.to_datetime(gdf['date'])  # Ensure 'time' is datetime
        group_col = 'date'
    else:
        raise ValueError("Neither 'year' nor 'date' column exists in the GeoDataFrame.")

    # Loop through unique values of the chosen column.
    for group_value in gdf[group_col].unique():
        if group_col == 'date':
            # Filter by exact time slice if 'time' is used.
            group_gdf = gdf[gdf['date'] == group_value]
            title_value = group_value.strftime('%Y-%m-%d')  # Format date for title
        else:
            # Filter by year if 'year' is used.
            group_gdf = gdf[gdf['year'] == group_value]
            title_value = group_value  # Year value for title
        
        # Plotting.
        print(f"Year {group_value} range: {group_gdf[var_of_interest].min()} to {group_gdf[var_of_interest].max()}")

        fig, ax = plt.subplots(figsize=(8, 8))
        plot = group_gdf.plot(
            column=var_of_interest, 
            ax=ax,
            markersize=0.5,
            cmap="mako_r",
            vmin=min_val,
            vmax=max_val
        )

        # Get the position of the plot axes
        pos = ax.get_position()

        # Create the colorbar at the top with EXACTLY the same width as the plot
        cbar_ax = fig.add_axes([pos.x0, pos.y0 + pos.height + 0.02, pos.width, 0.03])

        # Create the colorbar
        mappable = plot.get_children()[0]
        cbar = fig.colorbar(mappable, cax=cbar_ax, orientation='horizontal')

        # Set only endpoints (or endpoints + midpoint) on the colorbar
        if min_val is not None and max_val is not None:
            # Option 1: Only show endpoints
            cbar.set_ticks([min_val, max_val])
            
            # Option 2: Show endpoints and midpoint (uncomment if you want this)
            # midpoint = (min_val + max_val) / 2
            # cbar.set_ticks([min_val, midpoint, max_val])

        # Set the colorbar label
        cbar.set_label(f"{var_dict[var_of_interest]}")

        # Make sure ticks and label are on top
        cbar_ax.xaxis.set_ticks_position('top')
        cbar_ax.xaxis.set_label_position('top')

        # Set the plot options.
        #ax.set_title(f"{title_value#}", pad=50)
        ax.set_xticks([])
        ax.set_yticks([])

        # Don't use tight_layout() as it can adjust the carefully positioned axes
        # Instead, use a fixed adjustment to account for the colorbar
        plt.subplots_adjust(top=0.85)
        
        # # Pull colorbar units from var_dict.
        # colorbar = plot.get_figure().axes[-1]
        # colorbar.set_ylabel(f"{var_dict[var_of_interest]}")

        # # Set the plot options.
        # ax.set_title(f"{title_value}")
        # ax.set_xticks([])
        # ax.set_yticks([])
        # plt.tight_layout()

        # Save the plot.
        output_file = os.path.join(output_folder, f"{provider_name}-{var_of_interest}_{title_value}.png")
        plt.savefig(output_file, dpi=300, transparent=True)
        print(f'Saved figure to: {output_file}')
        plt.close()

def calculate_avg_df(df, dataset_name, time_column):

    if time_column == 'year':
        df['time'] = pd.to_datetime(df['year'], format='%Y')
    elif time_column == 'date':
        df['time'] = pd.to_datetime(df['date'], format='%Y-%m-%d')

    if dataset_name == "planet.forest_carbon_diligence":
        df_avg = df.groupby('time')['aboveground_live_carbon_density'].mean() / 0.476
        df_low = df.groupby('time')['aboveground_live_carbon_density_uncertainty_lower_bound'].mean() / 0.476
        df_high = df.groupby('time')['aboveground_live_carbon_density_uncertainty_upper_bound'].mean() / 0.476

        combined_df = pd.concat([df_avg, df_low, df_high], axis=1)
        combined_df.columns = ['average', 'lower_bound', 'upper_bound']
        combined_df.reset_index(inplace=True)
    
    elif dataset_name == "planet.forest_carbon_monitoring":
        df_avg = df.groupby('time')['aboveground_live_carbon_density'].mean() / 0.476
        df_low = df.groupby('time')['aboveground_live_carbon_density_uncertainty_lower_bound'].mean() / 0.476
        df_high = df.groupby('time')['aboveground_live_carbon_density_uncertainty_upper_bound'].mean() / 0.476

        combined_df = pd.concat([df_avg, df_low, df_high], axis=1)
        combined_df.columns = ['average', 'lower_bound', 'upper_bound']
        combined_df.reset_index(inplace=True)

    elif dataset_name in ["chloris.aboveground_biomass_stock_and_change_30_m", 
                          "chloris.aboveground_biomass_stock_and_change_10_m"]:
        df_avg = df.groupby('time')['aboveground_biomass_stock'].mean()
        df_error = df.groupby('time')['aboveground_biomass_stock_standard_error'].mean()

        combined_df = pd.concat([df_avg, df_error], axis=1)
        combined_df.columns = ['average', 'standard_error']
        combined_df.reset_index(inplace=True)

    elif dataset_name in ["kanop.monitoring_25_m",
                          "kanop.monitoring_10_m"]:
        combined_df = df.groupby('time')['living_aboveground_biomass'].mean()
        combined_df = combined_df.reset_index(name='living_aboveground_biomass')

    elif dataset_name == "kanop.screening_25_m":
        combined_df = df.groupby('time')['living_aboveground_biomass'].mean()
        combined_df = combined_df.reset_index(name='living_aboveground_biomass')

    return combined_df



def plot_dataset_timeseries(avg_df, dataset_name, output_folder):
    """
    Generates an average timeseries plot for a specific dataset with uncertainty.
    """
    os.makedirs(output_folder, exist_ok=True)

    if dataset_name == "planet.forest_carbon_diligence":
        filtered_df = avg_df.dropna(subset=[
            'time', 
            'average', 
            'lower_bound', 
            'upper_bound'
        ])
    elif dataset_name == "planet.forest_carbon_monitoring":
        filtered_df = avg_df.dropna(subset=[
            'time', 
            'average', 
            'lower_bound', 
            'upper_bound'
        ])
    elif dataset_name in ["chloris.aboveground_biomass_stock_and_change_30_m",
                          "chloris.aboveground_biomass_stock_and_change_10_m"]:
        filtered_df = avg_df.dropna(subset=[
            'time', 
            'average', 
            'standard_error'
        ])
    elif dataset_name in ["kanop.monitoring_25_m",
                          "kanop.monitoring_10_m"]:
        filtered_df = avg_df.dropna(subset=[
            'time', 
            'living_aboveground_biomass'
        ])
    elif dataset_name == "kanop.screening_25_m":
        filtered_df = avg_df.dropna(subset=[
            'time', 
            'living_aboveground_biomass'
        ])
    else:
        raise ValueError("Provider must be one of ['Planet', 'Chloris', 'Kanop']")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if dataset_name == "planet.forest_carbon_diligence":
        # Planet AGB with uncertainty
        ax.plot(filtered_df['time'], filtered_df['average'], label='Planet Diligence AGB (adjusted)', marker='o', color='#3F4E96')
        ax.fill_between(
            filtered_df['time'], 
            filtered_df['lower_bound'], 
            filtered_df['upper_bound'], 
            color='#3F4E96', alpha=0.2, label='Planet Diligence 90% Prediction Interval'
        )
    elif dataset_name == "planet.forest_carbon_monitoring":
        # Planet AGB with uncertainty
        ax.plot(filtered_df['time'], filtered_df['average'], label='Planet Monitoring AGB (adjusted)', marker='o', color='#E65832')
        ax.fill_between(
            filtered_df['time'], 
            filtered_df['lower_bound'], 
            filtered_df['upper_bound'], 
            color='#E65832', alpha=0.2, label='Planet Monitoring 90% Prediction Interval'
        )
    elif dataset_name in ["chloris.aboveground_biomass_stock_and_change_30_m",
                          "chloris.aboveground_biomass_stock_and_change_10_m"]:
        # Chloris AGB with uncertainty
        ax.plot(filtered_df['time'], filtered_df['average'], label='Chloris AGB', marker='o', color='#72924A')
        ax.fill_between(
            filtered_df['time'], 
            filtered_df['average'] - filtered_df['standard_error'], 
            filtered_df['average'] + filtered_df['standard_error'], 
            color='#72924A', alpha=0.2, label='Chloris Standard Error'
        )
    elif dataset_name in ["kanop.monitoring_25_m",
                          "kanop.monitoring_10_m",
                          "kanop.screening_25_m"]:
        # Kanop AGB (no uncertainty)
        ax.plot(filtered_df['time'], filtered_df['living_aboveground_biomass'], label='Kanop AGB', marker='o', color='#F09E3A')
    
    # Add labels, legend, and title
    ax.set_xlabel('Date')
    ax.set_ylabel('AGB (Mg/ha)')
    ax.legend(loc='upper right')
    ax.grid(False)
    
    # Save the plot
    output_file = os.path.join(output_folder, f"{dataset_name}_AGB_timeseries.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_customer_analysis(reprojections_to_query, data_requests_to_query,
                           dataset_mapping, variable_dicts, year, 
                           vars_of_interest, time_columns, output_dir):
    """ 
    Plots the raster image analysis for a customer. 
    
    The function first queries the data for each reprojection ID in the list
    reprojections_to_query and pulls all the data into a dictionary where each 
    individual provider dataframe is linked to its dataset name. It then loops
    through each dataset in the dictionary and plots 1) all of the variables for
    one year of interest (for the purposes of demonstrating each of the variables
    avaiable in each dataset), and 2) plots a full timeseries of a single variable
    of interest (e.g. for comparison with other data providers or examination of 
    change over time).

    Args:
    reprojections_to_query: a list of each reprojection ID to include in the analysis
    dataset_mapping: a dictionary linking dataset IDs with human readable names in the 
    format provider.dataset_name (e.g. 'chloris.aboveground_biomass_stock_and_change_30_m')
    variable_dicts: a dictionary of dictionaries (one dictionary for each dataset) with 
    all the variables and associated units.
    year: The year to use for plotting all the variables.
    vars_of_interest: a dictionary specifying the variable of interest for each of the 
    datasets being plotted (this is the variable that is plotted in the raster timeseries)

    Outputs:
    A series of plots saved to output_dir. 
    """

    data_dict = {}

    # Query data for all the reprojections listed.
    for reprojection_id in reprojections_to_query:
        df = query_all_data_analytics(reprojection_id, dataset_mapping)
        gdf = df_to_gdf(df, crs="EPSG:4326")
        dataset_name, aoi_name, aoi_id = get_reprojection_details_by_id(reprojection_id, dataset_mapping)

        # Set the output folder for the csv of data.
        #csv_output_dir = f"{output_dir}/{aoi_name}"
        #os.makedirs(csv_output_dir, exist_ok=True)
        #df.to_csv(f'{csv_output_dir}/{dataset_name}.csv', sep='\t', encoding='utf-8', index=False, header=True)

        # Add gdf to dictionary with reprojection_id as key.
        data_dict[reprojection_id] = gdf

    # Query the raw_db data for all the data requests listed.
    for data_request_id in data_requests_to_query:
        df = query_all_data_raw(data_request_id, dataset_mapping)
        gdf = df_to_gdf(df, crs="EPSG:4326")
        dataset_name, aoi_name, aoi_id = get_data_request_details_by_id(data_request_id, dataset_mapping)

        data_dict[data_request_id] = gdf


    # Loop through each dataset (reprojection by reprojection) and plot the variables.
    for data_id in data_dict.keys():
        if data_id in reprojections_to_query:
            dataset_name, aoi_name, aoi_id = get_reprojection_details_by_id(data_id, dataset_mapping)
        elif data_id in data_requests_to_query:
            dataset_name, aoi_name, aoi_id = get_data_request_details_by_id(data_id, dataset_mapping)
        else:
            raise ValueError(f"data_id {data_id} is not in either reprojections_to_query or data_requests_to_query.")

        provider_name = dataset_name.split(".")[0]
        print(f'Plotting data for aoi: {aoi_name}, dataset: {provider_name} - {dataset_name}')

        # Subset to just one dataset at a time for plotting.
        gdf_plot = data_dict[data_id]

        # Extract the relevant dict of variables and units.
        if dataset_name in variable_dicts:
            var_dict = variable_dicts[dataset_name]
        else:
            print(f"Dataset {dataset_name} not found in variable dictionaries.")

        # # Plot all variables for just one year.
        # output_folder = f"{output_dir}/{aoi_name}/{dataset_name}/variables"
        # plot_all_vars(gdf_plot, var_dict, year, provider_name, output_folder)

        # Select the variable of interest based on dataset name and plot the variable timeseries.
        var_of_interest = vars_of_interest[dataset_name]
        print(f'For dataset {dataset_name}, variable of interest is : {var_of_interest}')
        output_folder = f"{output_dir}/{aoi_name}/{dataset_name}/timeseries"
        plot_all_years(gdf_plot, var_dict, var_of_interest, provider_name, output_folder)

        
        # # Get the name of the column with time information based on dataset name.
        # time_column = time_columns[dataset_name]
        # print(f"Time column for {dataset_name}: {time_column}")

        # # Calculate the average dataframe and plot timeseries
        # timeseries_output_folder = f"{output_dir}/{aoi_name}/{dataset_name}/average_timeseries"
        # os.makedirs(timeseries_output_folder, exist_ok=True)
        # avg_df = calculate_avg_df(gdf_plot, dataset_name, time_column)
        # avg_df.to_csv(f'{timeseries_output_folder}/average_timeseries_AGB.csv')
        # plot_dataset_timeseries(avg_df, dataset_name, timeseries_output_folder)

    return data_dict



if __name__ == "__main__":

    # Pull dictionaries from config file to link dataset IDs, dataset names, variables, and units.
    with open('config.yml', 'r') as file:
        config = yaml.safe_load(file)
    dataset_mapping = config['dataset_mapping']
    variable_dicts = config['variables_units']
    time_columns = config['time_format']

    # Pull customer specific parameters from yml file.
    with open('customer-ymls/orbify-newsletter.yml', 'r') as file:
        customer_config = yaml.safe_load(file)
    year = customer_config['year']
    vars_of_interest = customer_config['vars_of_interest']
    reprojections_to_query = customer_config['reprojections_to_query']
    data_requests_to_query = customer_config['data_requests_to_query']

    output_dir = "outputs/orbify-newsletter"

    # Plot the charts
    plot_customer_analysis(reprojections_to_query, 
                           data_requests_to_query,
                           dataset_mapping, 
                           variable_dicts, 
                           year, 
                           vars_of_interest,
                           time_columns,
                           output_dir)


    


    
