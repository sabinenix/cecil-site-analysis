import cecil
import json
import geopandas as gpd
import pandas as pd
from shapely.geometry import shape



def show_linked_reprojections(dataset_mapping):
    """For the current organisation, list all the reprojections with their 
    relevant datasets, resolutions, and the AOI names and IDs."""
    client = cecil.Client()

    # List all reprojections, data requests and AOIs in an organisation.
    reprojections = client.list_reprojections()
    data_requests = client.list_data_requests()
    aois = client.list_aois()

    # Convert objects to dictionaries
    reprojection_dicts = [vars(r) for r in reprojections]
    data_request_dicts = [vars(dr) for dr in data_requests]
    aoi_dicts = [vars(aoi) for aoi in aois]

    # Convert lists of dictionaries to DataFrames
    reprojection_df = pd.DataFrame(reprojection_dicts)
    data_request_df = pd.DataFrame(data_request_dicts)
    aoi_df = pd.DataFrame(aoi_dicts)

    # Merge reprojections with data requests
    merged_df = reprojection_df.merge(
        data_request_df, 
        left_on='data_request_id', 
        right_on='id', 
        suffixes=('_reproj', '_data_req')
    )

    # Add dataset names using the dataset_mapping dictionary
    merged_df['dataset_name'] = merged_df['dataset_id'].map(dataset_mapping)

    # Merge with AOI information
    merged_with_aois = merged_df.merge(
        aoi_df[['id', 'name']], 
        left_on='aoi_id', 
        right_on='id', 
        suffixes=('', '_aoi')
    )

    # Rename and select relevant columns
    result_df = merged_with_aois.rename(columns={'name': 'aoi_name'})
    result_df = result_df[['id_reproj', 'dataset_name', 'resolution', 'aoi_name', 'aoi_id']]

    # Display the resulting DataFrame
    return result_df


def get_reprojection_details_by_id(reprojection_id, dataset_mapping):
    """For a given reprojection ID, return the dataset name, aoi name and ID."""
    
    # Get full dataframe linking reprojections to datasets and AOIs.
    result_df = show_linked_reprojections(dataset_mapping)

    # Filter to get just the reprojection ID.
    filtered = result_df[result_df['id_reproj'] == reprojection_id]

    if filtered.empty:
        return {"error": "Reprojection ID not found"}

    # Extract details of interest for querying data.
    dataset_name = filtered['dataset_name'].iloc[0]
    print(dataset_name)
    aoi_name = filtered['aoi_name'].iloc[0]
    print(aoi_name)
    aoi_id = filtered['aoi_id'].iloc[0]
    print(aoi_id)

    return dataset_name, aoi_name, aoi_id


def query_all_data(reprojection_id, dataset_mapping):
    """Pull all data for a given reprojection_id and aoi_id."""

    client = cecil.Client()

    print(f'Attempting to query data for {reprojection_id}')

    # Get details for querying the data from the reprojection_id.
    dataset_name, aoi_name, aoi_id = get_reprojection_details_by_id(reprojection_id, 
                                                          dataset_mapping)
    print(f'Pulling {dataset_name} for AOI {aoi_name} and AOI id {aoi_id}.')
    
    # Full table scan to pull all data for a given reprojection_id.
    df = client.query(f'''
        SELECT * 
        FROM 
            {dataset_name}
        WHERE
            aoi_id = '{aoi_id}' AND reprojection_id = '{reprojection_id}'
        ''')
    return df


def df_to_gdf(df, crs="EPSG:4326"):
    """
    Convert a DataFrame (from Snowflake) to GeoDataFrame with 'geometry' column.
    
    Args: 
    df: DataFrame pulled from Snowflake with 'pixel boundary' column.
    crs: (String) CRS of the data. Default is 'EPSG:4326'. 
    
    Returns:
    gdf: GeoDataFrame created using the 'geometry' column.
    """
    
    df["geometry"] = df["pixel_boundary"].apply(
    lambda x: shape(json.loads(x)) if isinstance(x, str) else None
    )
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=crs)
    return gdf
