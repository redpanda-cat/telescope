import React, { useState, useEffect } from "react";

import Paper from "@mui/material/Paper";
import Table from "@mui/material/Table";
import TableBody from "@mui/material/TableBody";
import TableCell from "@mui/material/TableCell";
import TableHead from "@mui/material/TableHead";
import TableRow from "@mui/material/TableRow";
import * as d3 from "d3";

const TableComponent = ({ filters }) => {
  const [genes, setGenes] = useState([]);

  useEffect(() => {
    if (filters) {
      const { timepoints, celltypes, clones } = filters;

      const HIGHLIGHT_URL = `/api/degenes?timepoint=${timepoints.join(
        ","
      )}&phenotype=${celltypes.join(",")}`;

      const URL = clones
        ? `${HIGHLIGHT_URL}&clone=${clones.join(",")}`
        : HIGHLIGHT_URL;

      fetch(URL)
        .then((res) => res.json())
        .then((result) => setGenes(result));
    } else {
      setGenes([]);
    }
  }, [filters]);

  return (
    <Paper sx={{ p: 2 }}>
      {filters ? (
        <GeneTable genes={genes} />
      ) : (
        "Please interact with the Sankey to see differentiated genes!"
      )}
    </Paper>
  );
};

const GeneTable = ({ genes }) => {
  return (
    <div style={{ height: 300, overflowY: "scroll", overflowX: "none" }}>
      <Table size="small">
        <TableHead>
          <TableRow>
            <TableCell>Gene</TableCell>
            <TableCell>Adjusted P value</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {genes.map((gene) => (
            <TableRow key={gene["gene"]}>
              <TableCell>{gene["gene"]}</TableCell>
              <TableCell>{d3.format("0.3")(gene["p"])}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
};

export default TableComponent;
