import React from "react";
import Paper from "@mui/material/Paper";
import Typography from "@mui/material/Typography";

const Info = ({ data, cells, filters }) => {
  let text = "None";
  if (filters) {
    const { timepoints, celltypes, clones } = filters;

    if (clones) {
      text = "Link : ";
      text += `${timepoints[0]} ${celltypes[0]} - ${timepoints[1]} ${celltypes[1]}`;
    } else {
      text = "Node : ";
      text += `${timepoints[0]} ${celltypes[0]}`;
    }
  }
  return (
    <Paper sx={{ p: 2 }}>
      <Typography variant="h5">
        Cells : {cells ? cells.length : data.length} / {data.length}
      </Typography>
      <Typography variant="h5">Selected</Typography>
      <Typography>{text}</Typography>
      {filters && filters.hasOwnProperty("clones") ? (
        <Typography>Clones : {filters.clones.join(", ")}</Typography>
      ) : null}
      <p></p>
    </Paper>
  );
};

export default Info;
