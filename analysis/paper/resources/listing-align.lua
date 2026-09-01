-- Left-align captioned code listings in docx output.
--
-- Quarto wraps a crossref-able float in a one-cell table so the caption stays
-- with its content, and gives that table an AlignCenter colspec. Pandoc's docx
-- writer turns cell alignment into direct <w:jc> on every paragraph inside,
-- which centres the code. Direct formatting beats paragraph styles, so this
-- cannot be corrected in the reference doc -- the colspec has to change.
--
-- Only wrappers holding a code block and nothing else are touched, so figure
-- and table floats stay centred.

if not FORMAT:match("docx") then
  return {}
end

local function listing_only(blocks)
  local code, other = false, false
  pandoc.walk_block(pandoc.Div(blocks), {
    CodeBlock = function() code = true end,
    Image = function() other = true end,
    Table = function() other = true end,
  })
  return code and not other
end

function Table(tbl)
  local body = tbl.bodies[1]
  if not body or not body.body[1] then return nil end
  local cell = body.body[1].cells[1]
  if not cell or not listing_only(cell.contents) then return nil end

  for i, spec in ipairs(tbl.colspecs) do
    tbl.colspecs[i] = { pandoc.AlignLeft, spec[2] }
  end
  return tbl
end
