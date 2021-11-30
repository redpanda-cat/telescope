import data from "./data.json";
import { Sankey } from "@shahlab/planetarium";

function App() {
  return (
    <Sankey
      data={data}
      width={800}
      height={700}
      subsetParam="cell_type"
      cloneParam="clone_id"
      timepointOrder={["Pre", "Post"]}
      timepointParam="treatment"
    />
  );
}

export default App;
