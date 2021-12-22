import React, { useState, useEffect } from "react";

import _ from "lodash";
import { scaleOrdinal, schemeTableau10 } from "d3";

import Box from "@mui/material/Box";
import Paper from "@mui/material/Paper";
import Grid from "@mui/material/Grid";
import { Sankey, UMAP, sortAlphanumeric } from "@shahlab/planetarium";
import Info from "./sankey/Info";
import Table from "./sankey/Table";

const App = ({ subsetParam, cloneParam, timepointParam, ...args }) => {
  const [data, setData] = useState(null);
  const [filters, setFilters] = useState(null);
  const [highlightedIDs, setHighlightedIDs] = useState(null);

  // Initial data fetch
  useEffect(() => {
    fetch("/api/data")
      .then((res) => res.json())
      .then((result) => {
        setData(result);
      });
  }, []);

  // Data fetch for highlighted cells
  useEffect(() => {
    if (filters) {
      const { timepoints, celltypes, clones } = filters;

      const HIGHLIGHT_URL = `/api/cells?timepoint=${timepoints.join(
        ","
      )}&phenotype=${celltypes.join(",")}`;

      const URL = clones
        ? `${HIGHLIGHT_URL}&clone=${clones.join(",")}`
        : HIGHLIGHT_URL;

      fetch(URL)
        .then((res) => res.json())
        .then((result) => setHighlightedIDs(result));
    } else {
      setHighlightedIDs(null);
    }
  }, [filters]);

  const onNodeClick = (node) => {
    if (!node) {
      setFilters(null);
    } else {
      const timepoints = [node[timepointParam]];
      const celltypes = [node[subsetParam]];
      setFilters({ timepoints, celltypes });
    }
  };

  const onLinkClick = (link) => {
    if (!link) {
      setFilters(null);
    } else {
      const { node0, node1 } = link;

      const timepoints = [node0[timepointParam], node1[timepointParam]];
      const celltypes = [node0[subsetParam], node1[subsetParam]];
      const clones = link[cloneParam];
      setFilters({ timepoints, celltypes, clones });
    }
  };

  const subsets = data
    ? _.uniq(data.map((datum) => datum[subsetParam])).sort(sortAlphanumeric)
    : [];
  const colors = scaleOrdinal(subsets, schemeTableau10);

  return data ? (
    <Box sx={{ p: 1 }}>
      <Grid
        container
        spacing={2}
        direction="column"
        sx={{ width: args.width * 2 + 75 }}
      >
        <Grid container spacing={2} item direction="row">
          <Grid item>
            <Paper>
              <Sankey
                data={data}
                subsetParam={subsetParam}
                cloneParam={cloneParam}
                timepointParam={timepointParam}
                onNodeClick={onNodeClick}
                onLinkClick={onLinkClick}
                colorScale={colors}
                {...args}
              />
            </Paper>
          </Grid>
          <Grid item>
            <Paper>
              <UMAP
                height={args.height}
                width={args.width}
                data={data}
                highlightIDs={highlightedIDs}
                xParam={"UMAP_1"}
                yParam={"UMAP_2"}
                subsetParam={subsetParam}
                idParam={"cell_id"}
                disable={true}
                colorScale={colors}
              />
            </Paper>
          </Grid>
        </Grid>
        <Grid container item direction="row" spacing={2}>
          <Grid item xs={4}>
            <Info data={data} cells={highlightedIDs} filters={filters} />
          </Grid>
          <Grid item xs={8}>
            <Paper>
              <Table filters={filters} />
            </Paper>
          </Grid>
        </Grid>
      </Grid>
      <div style={{ position: "absolute", bottom: 10, right: 10 }}>V3.0</div>
    </Box>
  ) : null;
};

const Dev = () => {
  return (
    <App
      width={800}
      height={700}
      subsetParam="cell_type"
      cloneParam="clone_id"
      timepointOrder={["Pre", "Post"]}
      timepointParam="treatment"
    />
  );
};

const Prod = () => {
  return <App {...window.appArgs} />;
};

export default process.env.NODE_ENV === "development" ? Dev : Prod;
