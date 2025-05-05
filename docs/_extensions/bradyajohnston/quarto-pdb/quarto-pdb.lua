
-- url-shortcode.lua
-- A Pandoc Lua filter that creates a shortcode for generating URLs

-- Function to create a URL from the given argument
function pdb_link(args, kwargs)
  local code = pandoc.utils.stringify(args[1])
  
  -- Create a link element with the URL as both text and destination
  return pandoc.Link(code, "https://rcsb.org/structure/" .. code)
end

-- Register the shortcode
return {
  ['pdb'] = pdb_link
}