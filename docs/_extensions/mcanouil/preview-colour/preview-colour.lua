--[[
# MIT License
#
# Copyright (c) Mickaël Canouil
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
]]

local function is_empty(s)
  return s == nil or s == ''
end

local preview_colour_meta = {
  ["text"] = true,
  ["code"] = true
}

function get_colour_preview_meta(meta)
  local preview_colour_text = true
  local preview_colour_code = true
  if not is_empty(meta['preview-colour']) then
    if not is_empty(meta['preview-colour']['text']) then
      preview_colour_text = meta['preview-colour']['text']
    end
    if not is_empty(meta['preview-colour']['code']) then
      preview_colour_code = meta['preview-colour']['code']
    end
  end
  meta['preview-colour'] = {
    ["text"] = preview_colour_text,
    ["code"] = preview_colour_code
  }
  preview_colour_meta = meta['preview-colour']
  return meta
end


function get_colour(element)
  function get_hex_color(n)
    return '#' .. string.rep('[0-9a-fA-F]', n)
  end

  local hex = nil
  for i = 6, 3, -1 do
    hex = element.text:match('(' .. get_hex_color(i) .. ')')
    if (i == 5 or i == 4) and hex ~= nil then
      hex = nil
      break
    end
    if hex ~= nil and (i == 6 or i == 3) then
      break
    end
  end
  if hex == nil then
    hex = element.text:match('(rgb%s*%(%s*%d+%s*,%s*%d+%s*,%s*%d+%s*%))')
  end
  if hex == nil then
    hex = element.text:match('(hsl%s*%(%s*%d+%s*,%s*%d+%s*%%,%s*%d+%s*%%s*%))')
  end

  return hex
end

function process_str(element, meta)
  if preview_colour_meta['text'] == false then
    return element
  end
  if preview_colour_meta['text'] == true then
    hex = get_colour(element)
    if hex ~= nil then
      if quarto.doc.is_format("html:js") then
        colour_preview_mark = "<span style=\"display: inline-block; color: " .. hex .. ";\">&#9673;</span>"
        new_text = string.gsub(
          element.text,
          hex,
          hex .. colour_preview_mark
        )
        return pandoc.RawInline('html', new_text)
      elseif quarto.doc.is_format("latex") then
        colour_preview_mark = "\\textcolor[HTML]{" .. string.gsub(hex, '#', '') .. "}{\\textbullet}"
        new_text = string.gsub(
          element.text,
          hex,
          "\\" .. hex .. colour_preview_mark
        )
        return pandoc.RawInline('latex', new_text)
      end
    end
  end
end

function process_code(element)
  if preview_colour_meta['code'] == false then
    return element
  end
  if preview_colour_meta['code'] == true then
    hex = get_colour(element)
    if hex ~= nil then
      if quarto.doc.is_format("html:js") then
        colour_preview_mark = "<span style=\"display: inline-block; color: " .. hex .. ";\">&#9673;</span>"
        return pandoc.Span({element, pandoc.RawInline('html', colour_preview_mark)})
      elseif quarto.doc.is_format("latex") then
        colour_preview_mark = "\\textcolor[HTML]{" .. string.gsub(hex, '#', '') .. "}{\\textbullet}"
        return pandoc.Span({element, pandoc.RawInline('latex', colour_preview_mark)})
      end
    end
  end
end

return {
  {Meta = get_colour_preview_meta},
  {Str = process_str},
  {Code = process_code}
}
