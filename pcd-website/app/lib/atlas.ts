import { PCDAtlas } from "../types";

export const ATLAS_RAW_URL =
  "https://raw.githubusercontent.com/aesenthilvanan-coder/pcd-atlas-data/main/PCD_global_atlas.json";

export async function fetchAtlas(): Promise<PCDAtlas> {
  const res = await fetch(ATLAS_RAW_URL, {
    cache: "no-store",          // always fresh — never stale cache
    next: { revalidate: 0 },    // Next.js: do not cache
  });
  if (!res.ok) throw new Error(`Atlas fetch failed: ${res.status}`);
  return res.json();
}
