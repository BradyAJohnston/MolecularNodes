function ensureLatexDeps()
  quarto.doc.use_latex_package("academicons")
end

local function ensureHtmlDeps()
  quarto.doc.add_html_dependency({
    name = "academicons",
    version = "1.9.4",
    stylesheets = { "assets/css/all.css", "assets/css/size.css" }
  })
end

local function isEmpty(s)
  return s == nil or s == ''
end

local function isValidSize(size)
  local validSizes = {
    "tiny",
    "scriptsize",
    "footnotesize",
    "small",
    "normalsize",
    "large",
    "Large",
    "LARGE",
    "huge",
    "Huge"
  }
  for _, v in ipairs(validSizes) do
    if v == size then
      return size
    end
  end
  return ""
end

local function convertToLatexSize(size)
  local sizeMap = {
    ["2xs"] = "tiny",
    ["xs"] = "scriptsize",
    ["sm"] = "small",
    ["lg"] = "large",
    ["xl"] = "Large",
    ["2xl"] = "huge",
    ["1x"] = "normalsize",
    ["2x"] = "huge",
    ["3x"] = "Huge"
  }
  return sizeMap[size] or size
end
local function convertToHtmlSize(size)
  local sizeMap = {
    ["tiny"] = "2xs",
    ["scriptsize"] = "xs",
    ["footnotesize"] = "xs",
    ["small"] = "sm",
    ["normalsize"] = "1x",
    ["large"] = "lg",
    ["Large"] = "xl",
    ["LARGE"] = "2xl",
    ["huge"] = "2xl",
    ["Huge"] = "3x"
  }
  return sizeMap[size] or size
end

return {
  ["ai"] = function(args, kwargs)
    local group = ""
    local icon = pandoc.utils.stringify(args[1])
    if #args > 1 then
      group = icon
      icon = pandoc.utils.stringify(args[2])
    end

    local color = pandoc.utils.stringify(kwargs["color"])
    local hcolor = pandoc.utils.stringify(kwargs["hcolor"])
    local pcolor = pandoc.utils.stringify(kwargs["pcolor"])

    local label = pandoc.utils.stringify(kwargs["label"])
    if isEmpty(label) then
      label = " aria-label=\"" .. icon .. "\""
    else
      label = " aria-label=\"" .. label .. "\""
    end

    local title = pandoc.utils.stringify(kwargs["title"])
    if not isEmpty(title) then
      title = " title=\"" .. title .. "\""
    end

    local size = pandoc.utils.stringify(kwargs["size"])
    local hsize = pandoc.utils.stringify(kwargs["hsize"])
    local psize = pandoc.utils.stringify(kwargs["psize"])


    -- detect html (excluding epub)
    if quarto.doc.isFormat("html:js") then
      ensureHtmlDeps()
      if not isEmpty(hsize) then
        size = hsize
      end
      size = convertToHtmlSize(size)
      if not isEmpty(size) then
        size = " ai-" .. size
      end
      if not isEmpty(hcolor) then
        color = hcolor
      end
      if not isEmpty(color) then
        color = " style=\"color:" .. color .. "\""
      end
      return pandoc.RawInline(
        'html',
        "<i class=\"ai " .. group .. " ai-" .. icon .. size .. "\"" .. title .. color .. label .. "></i>"
      )
      -- detect pdf / beamer / latex / etc
    elseif quarto.doc.is_format("pdf") then
      ensureLatexDeps()
      if not isEmpty(psize) then
        size = psize
      end
      size = isValidSize(convertToLatexSize(size))
      if not isEmpty(size) then
        size = "\\" .. size
      end
      if not isEmpty(pcolor) then
        color = pcolor
      end
      if not isEmpty(color) then
        color = "\\color{" .. color .. "}"
      end
      icon = icon:gsub("-square", "_SQUARE"):gsub("-", ""):gsub("_SQUARE", "-square")

      if isEmpty(size) and isEmpty(color) then
        return pandoc.RawInline('tex', "\\aiicon{" .. icon .. "}")
      else
        return pandoc.RawInline('tex', "{" .. size .. color .. "\\aiicon{" .. icon .. "}}")
      end
    else
      return pandoc.Null()
    end
  end
}
