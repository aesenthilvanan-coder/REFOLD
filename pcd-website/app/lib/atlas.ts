import { PCDAtlas } from "../types";

const ATLAS_URL =
  "https://raw.githubusercontent.com/aesenthilvanan-coder/pcd-atlas-data/main/PCD_global_atlas.json";

export async function fetchAtlas(): Promise<PCDAtlas> {
  const res = await fetch(ATLAS_URL, {
    next: { revalidate: 300 }, // cache for 5 min; daemon pushes ~every 45s
  });
  if (!res.ok) throw new Error(`Atlas fetch failed: ${res.status}`);
  return res.json();
}
