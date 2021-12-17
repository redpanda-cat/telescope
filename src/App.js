import devData from "./data.json";
import { Sankey } from "@shahlab/planetarium";

const App = (args) => {
  return (
    <div>
      <Sankey {...args} />
      <div style={{ position: "absolute", bottom: 10, right: 10 }}>V2.0</div>
    </div>
  );
};

const Dev = () => {
  return (
    <App
      data={devData}
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
