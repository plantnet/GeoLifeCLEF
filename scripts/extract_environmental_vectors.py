import argparse
from pathlib import Path

import pandas as pd
import numpy as np
from tqdm import tqdm

from data_loading.environmental_raster import PatchExtractor


def compute_environmental_vectors(df, extractor, as_dataframe=True):
    def compute_environmental_vector(index):
        position = df.iloc[index][["latitude", "longitude"]]
        return extractor[position]

    environmental_vectors = [compute_environmental_vector(i) for i in tqdm(range(len(df)))]
    environmental_vectors = np.asarray(environmental_vectors)

    if as_dataframe:
        raster_names = [raster.name for raster in extractor.rasters_fr]
        df_env = pd.DataFrame(environmental_vectors, columns=raster_names)
        df_env.index = df.index
        df_env.index.name = "observation_id"
        return df_env
    else:
        return environmental_vectors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts the environmental vectors from the rasters"
    )
    parser.add_argument(
        "data_path",
        type=Path,
        help="path to directory containing the data",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="output CSV file containing the environmental vectors",
    )
    args = parser.parse_args()

    # Loading observations
    df_fr_train = pd.read_csv(args.data_path / "observations" / "observations_fr_train.csv", sep=";", index_col="observation_id")
    df_us_train = pd.read_csv(args.data_path / "observations" / "observations_us_train.csv", sep=";", index_col="observation_id")
    df_train = pd.concat((df_fr_train, df_us_train))

    df_fr_test = pd.read_csv(args.data_path / "observations" / "observations_fr_test.csv", sep=";", index_col="observation_id")
    df_us_test = pd.read_csv(args.data_path / "observations" / "observations_us_test.csv", sep=";", index_col="observation_id")
    df_test = pd.concat((df_fr_test, df_us_test))

    df_train = df_train[["latitude", "longitude"]]
    df = pd.concat((df_train, df_test))

    # Loading rasters
    print("Loading rasters...")
    extractor = PatchExtractor(args.data_path / "rasters", size=1)
    extractor.add_all_rasters(nan=np.nan)

    # Computing the environmental vectors
    print("Extracting environmental vectors...")
    df_env = compute_environmental_vectors(df, extractor, as_dataframe=True)
    df_env.sort_index(inplace=True)

    # Saving the result to CSV file
    df_env.to_csv(args.output_file, sep=";")
