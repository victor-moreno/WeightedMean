#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom jsonlite fromJSON toJSON
NULL

weightedMeanClass <- R6::R6Class(
    "weightedMeanClass",
    inherit = weightedMeanBase,
    private = list(

        .run = function() {

            # ── 0. Guard: need at least one variable ──────────────────────────
            if (is.null(self$options$vars) || length(self$options$vars) == 0) {
                self$results$propTableWidget$setContent(
                    '<div style="padding:16px;color:#888;font-family:sans-serif;font-size:13px;">
                     Select variables and (optionally) provide a property CSV path to begin.</div>'
                )
                return()
            }

            selectedVars <- self$options$vars

            # ── 1. Read the property CSV ──────────────────────────────────────
            propPath     <- self$options$propFile
            propDF       <- NULL
            fileWarnHtml <- ""

            if (!is.null(propPath) && nchar(trimws(propPath)) > 0) {
                propPath <- path.expand(trimws(propPath))
                if (file.exists(propPath)) {
                    tryCatch({
                        propDF <- read.csv(propPath, stringsAsFactors = FALSE)
                        names(propDF) <- trimws(names(propDF))
                        required <- c("Var", "weightX", "weightY", "weightZ")
                        missing  <- setdiff(required, names(propDF))
                        if (length(missing) > 0)
                            jmvcore::reject(paste("CSV missing columns:",
                                                  paste(missing, collapse = ", ")))
                    }, error = function(e) {
                        jmvcore::reject(paste("Could not read CSV:", e$message))
                    })
                } else {
                    fileWarnHtml <- paste0(
                        '<div style="padding:8px 10px;color:#c62828;font-family:sans-serif;',
                        'font-size:12px;border:1px solid #ef9a9a;border-radius:4px;',
                        'background:#ffebee;margin-bottom:8px;">',
                        '&#9888; File not found: <code>', propPath,
                        '</code> &mdash; using default weights (1).</div>'
                    )
                }
            }

            # No CSV / not found: skeleton with weight = 1
            if (is.null(propDF)) {
                propDF <- data.frame(
                    Var     = selectedVars,
                    weightX = rep(1, length(selectedVars)),
                    weightY = rep(1, length(selectedVars)),
                    weightZ = rep(1, length(selectedVars)),
                    stringsAsFactors = FALSE
                )
            }

            # ── 2. Apply user edits from hidden JSON option ───────────────────
            overridesJSON <- self$options$tableOverrides
            # Guard: treat empty / whitespace-only / "{}" / NULL as no overrides
            hasOvr <- !is.null(overridesJSON) &&
                      nchar(trimws(overridesJSON)) > 2 &&
                      trimws(overridesJSON) != "{}"

            if (hasOvr) {
                tryCatch({
                    overrides <- jsonlite::fromJSON(overridesJSON,
                                                   simplifyDataFrame = FALSE)
                    for (varName in names(overrides)) {
                        idx <- which(propDF$Var == varName)
                        if (length(idx) == 0) next
                        varOvr <- overrides[[varName]]
                        for (col in names(varOvr)) {
                            if (col %in% names(propDF))
                                propDF[idx, col] <- as.numeric(varOvr[[col]])
                        }
                    }
                }, error = function(e) { })
            }

            # ── 3. Match table ────────────────────────────────────────────────
            matchTable <- self$results$matchInfo
            matchTable$deleteRows()
            for (v in selectedVars) {
                idx <- which(propDF$Var == v)
                if (length(idx) > 0) {
                    matchTable$addRow(rowKey = v, values = list(
                        variable = v, matched = "\u2713 Yes",
                        weightX  = propDF$weightX[idx[1]],
                        weightY  = propDF$weightY[idx[1]],
                        weightZ  = propDF$weightZ[idx[1]]
                    ))
                } else {
                    matchTable$addRow(rowKey = v, values = list(
                        variable = v, matched = "\u2717 No (weight=1)",
                        weightX  = 1, weightY = 1, weightZ = 1
                    ))
                }
            }

            # ── 4. Build editable HTML grid ───────────────────────────────────
            displayDF <- propDF[propDF$Var %in% selectedVars, , drop = FALSE]
            unmatched <- setdiff(selectedVars, propDF$Var)
            if (length(unmatched) > 0) {
                displayDF <- rbind(displayDF, data.frame(
                    Var = unmatched, weightX = 1, weightY = 1, weightZ = 1,
                    stringsAsFactors = FALSE
                ))
            }
            ord <- match(selectedVars, displayDF$Var)
            displayDF <- displayDF[ord[!is.na(ord)], , drop = FALSE]

            self$results$propTableWidget$setContent(
                paste0(fileWarnHtml,
                       private$.buildEditableGrid(displayDF, overridesJSON))
            )

            # ── 5. Compute weighted means ─────────────────────────────────────
            dataDF    <- self$data
            weightCol <- self$options$weightDim

            resultsTable <- self$results$results
            resultsTable$deleteRows()
            weightedValues <- numeric(0)
            weightValues   <- numeric(0)

            for (v in selectedVars) {
                if (!v %in% names(dataDF)) next
                col <- as.numeric(dataDF[[v]])
                col <- col[!is.na(col)]
                mn  <- mean(col)
                idx <- which(propDF$Var == v)
                w   <- as.numeric(if (length(idx) > 0) propDF[idx[1], weightCol] else 1)
                resultsTable$addRow(rowKey = v, values = list(
                    variable     = v, mean = mn,
                    weightUsed   = w, weightedMean = mn * w
                ))
                weightedValues <- c(weightedValues, mn)
                weightValues   <- c(weightValues, w)
            }

            # ── 6. Overall weighted mean ──────────────────────────────────────
            overallTable <- self$results$overallResult
            overallTable$deleteRows()
            if (length(weightedValues) > 0 && sum(weightValues) != 0) {
                overallTable$addRow(rowKey = "overall", values = list(
                    dimension           = weightCol,
                    overallWeightedMean = sum(weightedValues * weightValues) /
                                         sum(weightValues)
                ))
            }
        },  # end .run

        # ── Helper: zero-dependency contenteditable HTML table ────────────────
        #
        # JS → R bridge strategy:
        #
        #   The postMessage/window.jamovi approach is unreliable because the
        #   exact message format and API availability varies across jamovi
        #   versions and is not publicly documented. We try every known variant,
        #   BUT we always show a copy-paste fallback so the user is never stuck:
        #
        #   1. Try all postMessage formats and window.jamovi APIs (auto)
        #   2. Always render the JSON in a read-only textarea after Save
        #   3. User pastes that JSON into the "Table Overrides (JSON)" TextBox
        #      in the Options panel — R reads it on next run
        #
        #   The tableOverrides option is therefore made VISIBLE in .u.yaml
        #   (not hidden) so the user can paste into it directly.
        #
        .buildEditableGrid = function(df, currentOverridesJSON) {

            nRows <- nrow(df)
            uid   <- paste0(sample(c(letters, 0:9), 8, replace = TRUE), collapse = "")

            # Build rows in R — values baked in, no JS data bootstrap
            rowsHtml <- paste(vapply(seq_len(nRows), function(i) {
                v  <- htmlEscape_(df$Var[i])
                wx <- formatC(as.numeric(df$weightX[i]), format = "f", digits = 4)
                wy <- formatC(as.numeric(df$weightY[i]), format = "f", digits = 4)
                wz <- formatC(as.numeric(df$weightZ[i]), format = "f", digits = 4)
                paste0(
                    '<tr data-var="', v, '">',
                    '<td class="wm-varname">', v, '</td>',
                    '<td class="wm-w" contenteditable="true">', wx, '</td>',
                    '<td class="wm-w" contenteditable="true">', wy, '</td>',
                    '<td class="wm-w" contenteditable="true">', wz, '</td>',
                    '</tr>'
                )
            }, character(1)), collapse = "\n")

            paste0(
'<style>
.wm-wrap { font-family:"Segoe UI",Arial,sans-serif; padding:2px 0 14px 0; }
.wm-wrap h4 { margin:0 0 10px 0; font-size:13px; color:#333; font-weight:600; }
.wm-tbl {
  border-collapse:collapse; font-size:12px; width:auto; min-width:340px;
}
.wm-tbl th {
  background:#e8edf5; padding:5px 14px; text-align:center;
  font-weight:600; font-size:11px; border:1px solid #c5cdd8;
  letter-spacing:.3px; white-space:nowrap;
}
.wm-tbl th:first-child { text-align:left; }
.wm-varname {
  padding:4px 14px; border:1px solid #dde1e7;
  background:#f5f6f8; color:#444;
  font-family:monospace; font-size:12px;
  white-space:nowrap; user-select:none;
}
.wm-w {
  padding:4px 10px; border:1px solid #dde1e7;
  text-align:right; font-family:monospace; font-size:12px;
  min-width:80px; outline:none; background:#fff;
}
.wm-w:focus {
  background:#fffde7; border-color:#1976d2;
  box-shadow:inset 0 0 0 1px #1976d2;
}
.wm-w.wm-dirty { background:#fff8e1; }
#wm-save-', uid, ' {
  margin-top:10px; padding:5px 16px; font-size:12px; cursor:pointer;
  background:#1976d2; color:#fff; border:none; border-radius:4px;
  font-family:"Segoe UI",Arial,sans-serif;
}
#wm-save-', uid, ':hover     { background:#1565c0; }
#wm-save-', uid, '.wm-unsaved { background:#e65100; }
#wm-bar-', uid, ' {
  display:none; margin-top:7px; padding:5px 10px; border-radius:4px;
  font-size:12px; background:#e8f5e9; color:#2e7d32; border:1px solid #a5d6a7;
}
.wm-note { font-size:11px; color:#999; margin-top:6px; }
/* ── Copy-paste fallback panel ── */
#wm-fallback-', uid, ' {
  display:none; margin-top:14px; padding:10px 12px;
  border:1px solid #bbdefb; border-radius:4px; background:#e3f2fd;
  font-size:12px; font-family:"Segoe UI",Arial,sans-serif;
}
#wm-fallback-', uid, ' strong { color:#0d47a1; }
#wm-ta-', uid, ' {
  display:block; width:100%; margin-top:6px; padding:5px 7px;
  font-family:monospace; font-size:11px;
  border:1px solid #90caf9; border-radius:3px;
  background:#fff; resize:vertical; min-height:46px;
  box-sizing:border-box;
}
#wm-cbtn-', uid, ' {
  margin-top:6px; padding:3px 10px; font-size:11px; cursor:pointer;
  background:#1565c0; color:#fff; border:none; border-radius:3px;
  font-family:"Segoe UI",Arial,sans-serif;
}
</style>

<div class="wm-wrap">
  <h4>Variable Properties &#8212; click a weight cell to edit</h4>
  <table class="wm-tbl" id="wm-tbl-', uid, '">
    <thead>
      <tr>
        <th>Variable</th>
        <th>weightX</th>
        <th>weightY</th>
        <th>weightZ</th>
      </tr>
    </thead>
    <tbody>
', rowsHtml, '
    </tbody>
  </table>

  <button id="wm-save-', uid, '" onclick="wmSave_', uid, '()">
    &#128190; Save Changes
  </button>
  <div id="wm-bar-', uid, '">
    &#10003; Changes sent &#8212; if results do not update, use the fallback below.
  </div>

  <!-- Copy-paste fallback — always shown after Save so user is never stuck -->
  <div id="wm-fallback-', uid, '">
    <strong>Apply weights:</strong>
    Copy the JSON below and paste it into the
    <em>Table Overrides (JSON)</em> field in the Options panel,
    then click away to re-run.<br>
    <textarea id="wm-ta-', uid, '" readonly></textarea>
    <button id="wm-cbtn-', uid, '" onclick="wmCopy_', uid, '()">
      &#128203; Copy
    </button>
  </div>

  <p class="wm-note">
    Click a weight cell and type to change its value.
    Press Enter or Tab to confirm. Variable names are read-only.
  </p>
</div>

<script>
(function () {
  var tbl = document.getElementById("wm-tbl-', uid, '");
  var btn = document.getElementById("wm-save-', uid, '");
  var bar = document.getElementById("wm-bar-', uid, '");

  tbl.addEventListener("input", function (e) {
    if (!e.target.classList.contains("wm-w")) return;
    e.target.classList.add("wm-dirty");
    btn.classList.add("wm-unsaved");
    btn.textContent = "\u26A0 Unsaved \u2014 click to save";
  });

  tbl.addEventListener("keydown", function (e) {
    if (!e.target.classList.contains("wm-w")) return;
    if (e.key === "Enter") { e.preventDefault(); e.target.blur(); }
    if (e.key === "Tab") {
      e.preventDefault();
      var cells = Array.prototype.slice.call(tbl.querySelectorAll(".wm-w"));
      var next  = cells[cells.indexOf(e.target) + (e.shiftKey ? -1 : 1)];
      if (next) next.focus();
    }
  });

  window["wmSave_', uid, '"] = function () {
    var overrides = {};
    tbl.querySelectorAll("tbody tr").forEach(function (tr) {
      var v  = tr.getAttribute("data-var");
      var cs = tr.querySelectorAll(".wm-w");
      var wx = parseFloat(cs[0].textContent.trim());
      var wy = parseFloat(cs[1].textContent.trim());
      var wz = parseFloat(cs[2].textContent.trim());
      if (!v || isNaN(wx) || isNaN(wy) || isNaN(wz)) return;
      overrides[v] = { weightX: wx, weightY: wy, weightZ: wz };
    });
    var json = JSON.stringify(overrides);

    /* ── Attempt every known jamovi JS→option bridge ── */
    /* Format 1: documented in some jamovi versions */
    try { window.parent.postMessage(
            { type: "setOption", name: "tableOverrides", value: json }, "*");
    } catch(e) {}
    /* Format 2: alternative format */
    try { window.parent.postMessage(
            { action: "setOption", name: "tableOverrides", value: json }, "*");
    } catch(e) {}
    /* Format 3: jamovi direct API */
    try {
      if (window.jamovi) {
        if (typeof window.jamovi.setOption === "function")
          window.jamovi.setOption("tableOverrides", json);
        if (window.jamovi.results &&
            typeof window.jamovi.results.setOption === "function")
          window.jamovi.results.setOption("tableOverrides", json);
      }
    } catch(e) {}

    btn.classList.remove("wm-unsaved");
    btn.textContent = "\uD83D\uDCBE Save Changes";
    bar.style.display = "block";
    setTimeout(function () { bar.style.display = "none"; }, 4000);

    /* Always populate the fallback textarea */
    document.getElementById("wm-ta-', uid, '").value = json;
    document.getElementById("wm-fallback-', uid, '").style.display = "block";
  };

  window["wmCopy_', uid, '"] = function () {
    var ta = document.getElementById("wm-ta-', uid, '");
    ta.select();
    try { document.execCommand("copy"); } catch(e) {}
    var cb = document.getElementById("wm-cbtn-', uid, '");
    cb.textContent = "\u2713 Copied!";
    setTimeout(function () { cb.textContent = "\uD83D\uDCCB Copy"; }, 2000);
  };
}());
</script>'
            )
        }   # end .buildEditableGrid
    )       # end private
)

# ── Module-level helper ───────────────────────────────────────────────────────
htmlEscape_ <- function(x) {
    x <- gsub("&",  "&amp;",  x, fixed = TRUE)
    x <- gsub("<",  "&lt;",   x, fixed = TRUE)
    x <- gsub(">",  "&gt;",   x, fixed = TRUE)
    x <- gsub('"',  "&quot;", x, fixed = TRUE)
    x
}