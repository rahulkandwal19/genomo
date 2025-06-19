
document.addEventListener('DOMContentLoaded', function() {
  const analysisForm = document.getElementById("analysisForm");
  const resultsSection = document.getElementById("resultsSection");
  const analysisChartCtx = document.getElementById("analysisChart").getContext('2d');
  const summaryCard = document.getElementById("summaryCard");
  const matchesTable = document.getElementById("matchesTable").querySelector('tbody');
  const loadingSpinner = document.getElementById("loadingSpinner");
  const submitBtn = document.getElementById("submitBtn");
  const resetBtn = document.getElementById("resetBtn");

  // Initialize chart
  const analysisChart = new Chart(analysisChartCtx, {
    type: 'bar',
    data: { labels: [], datasets: [] },
    options: {
      responsive: true,
      scales: { y: { beginAtZero: true } }
    }
  });

  // Form submission handler
  analysisForm.addEventListener("submit", async function(e) {
    e.preventDefault();
    
    try {
      // Show loading state
      submitBtn.disabled = true;
      submitBtn.innerHTML = '<i class="fas fa-spinner fa-spin mr-2"></i> Analyzing...';
      loadingSpinner.classList.remove('hidden');
      resultsSection.classList.add('hidden');

      // Get form data
      const formData = new FormData(e.target);
      const data = Object.fromEntries(formData.entries());

      // Call your actual backend API
      const response = await fetch("http://localhost:18080/analyze", {
        method: "POST",
        headers: {
          "Content-Type": "application/json"
        },
        body: JSON.stringify(data)
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const result = await response.json();
      console.log(result);
      displayResults(result);

    } catch (error) {
      console.error("Analysis failed:", error);
      showError(error.message || "Analysis failed. Please try again.");
    } finally {
      // Reset loading state
      submitBtn.disabled = false;
      submitBtn.innerHTML = '<i class="fas fa-play mr-2"></i> Run Analysis';
      loadingSpinner.classList.add('hidden');
    }
  });

  // Reset button handler
  resetBtn.addEventListener("click", function() {
    analysisForm.reset();
    resultsSection.classList.add('hidden');
  });

// function to display pathogen paths
function displayPathogenPaths(paths) {
  const pathsTable = document.getElementById("pathsTable").querySelector('tbody');
  pathsTable.innerHTML = '';
  
  if (paths && paths.length > 0) {
    paths.forEach((path, index) => {
      const row = document.createElement('tr');
      row.className = 'hover:bg-gray-50 transition';
      
      // Create path display with arrows
      let pathDisplay = '';
      path.forEach((node, i) => {
        pathDisplay += `<span class="bg-gray-100 px-2 py-1 rounded-md font-medium">${node}</span>`;
        if (i < path.length - 1) {
          pathDisplay += ' <i class="fas fa-arrow-right text-gray-400 mx-1"></i> ';
        }
      });
      
      row.innerHTML = `
        <td class="px-6 py-4 whitespace-nowrap font-medium">${index + 1}</td>
        <td class="px-6 py-4">
          <div class="flex items-center flex-wrap">${pathDisplay}</div>
        </td>
        <td class="px-6 py-4 whitespace-nowrap">
          <span class="inline-flex items-center px-3 py-1 rounded-full text-xs font-medium bg-blue-100 text-blue-800">
            <i class="fas fa-arrows-alt-h mr-1"></i>
            ${path.length - 1} steps
          </span>
        </td>
      `;
      pathsTable.appendChild(row);
    });
  } else {
    pathsTable.innerHTML = `
      <tr>
        <td colspan="3" class="px-6 py-8 text-center">
          <div class="flex flex-col items-center text-gray-500">
            <i class="fas fa-ban text-3xl mb-2"></i>
            <p class="font-medium">No pathogen transmission paths found</p>
          </div>
        </td>
      </tr>
    `;
  }
}

function displayResults(result) {
    // Parse summary result to extract the values 
    const matches = [];
    let pathogenName = "Unknown";
    let pathogenData = "";
    
    // Extract pathogen name if available in the result
    if (result.pathogen_name) {
        pathogenName = result.pathogen_name;
    } else if (result.summary && result.summary.includes("|P1|")) {
        // Parse from GDN format
        const pathogenMatch = result.summary.match(/\d+\|(P\d+)\|/);
        if (pathogenMatch) {
            pathogenName = pathogenMatch[1];
        }
    }


    if (result.summary) {
        const matchLines = result.summary.split('\n');
        for (const line of matchLines) {
            const match = line.match(/•\s+(\S+)\s+→\s+Score:\s+([\d.]+)/);
            if (match) {
                matches.push({
                    position: match[1].trim(),
                    type: "CDS Match",
                    score: parseFloat(match[2]),
                    //pathogen: pathogenName,
                    sequence: "" 
                });
            }
        }
    }   
    // Sort matches by score (highest first)
    matches.sort((a, b) => b.score - a.score);

    // Update summary card
    summaryCard.innerHTML = `
        <h3 class="font-medium text-indigo-800 mb-3">Analysis Summary</h3>
        <div class="space-y-3">
            <div class="flex justify-between">
                <span class="text-gray-600">Pathogen Identified</span>
                <span class="font-medium flex items-center">
                    <i class="fas fa-virus text-red-500 mr-2"></i>
                    ${pathogenName}
                </span>
            </div>
            <div class="flex justify-between">
                <span class="text-gray-600">Total Matches</span>
                <span class="font-medium ${matches.length > 0 ? 'text-blue-600' : 'text-gray-500'}">
                    ${matches.length}
                </span>
            </div>
            <div class="flex justify-between">
                <span class="text-gray-600">Highest Score</span>
                <span class="font-medium ${getScoreColorClass(matches[0]?.score || 0)}">
                    ${matches.length > 0 ? matches[0].score.toFixed(3) : 'N/A'}
                </span>
            </div>
            <div class="flex justify-between">
                <span class="text-gray-600">Risk Assessment</span>
                <span class="font-medium ${getRiskLevelClass(matches)}">
                    ${determineRiskLevel(matches)}
                </span>
            </div>
        </div>
    `;

    // Update chart
    analysisChart.data.labels = matches.map(m => m.position);
    analysisChart.data.datasets = [{
        label: 'Match Confidence',
        data: matches.map(m => m.score),
        backgroundColor: matches.map(m => getScoreColor(m.score, 'background')),
        borderColor: matches.map(m => getScoreColor(m.score, 'border')),
        borderWidth: 2,
        borderRadius: 4
    }];
    analysisChart.options.scales.y = {
        min: 0,
        max: 1.0,
        ticks: {
            callback: function(value) {
                return value.toFixed(1);
            }
        }
    };
    analysisChart.update();

    // Update matches table
    matchesTable.innerHTML = '';
    if (matches.length > 0) {
        matches.forEach(match => {
            const row = document.createElement('tr');
            row.className = 'hover:bg-gray-50 transition';
            row.innerHTML = `
                <td class="px-6 py-4 whitespace-nowrap font-mono">
                    ${match.position}
                </td>
                <td class="px-6 py-4 whitespace-nowrap">
                    <span class="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium ${getScoreColorClass(match.score)}">
                        ${match.type}
                    </span>
                </td>
                <td class="px-6 py-4 whitespace-nowrap">
                    <div class="flex items-center">
                        <div class="h-2 w-full bg-gray-200 rounded-full overflow-hidden">
                            <div class="h-2 ${getScoreColor(match.score, 'bar')}" 
                                 style="width: ${match.score * 100}%"></div>
                        </div>
                        <span class="ml-2 font-mono text-sm ${getScoreTextClass(match.score)}">
                            ${match.score.toFixed(3)}
                        </span>
                    </div>
                </td>
             
            `;
            matchesTable.appendChild(row);
        });
    } else {
        matchesTable.innerHTML = `
            <tr>
                <td colspan="5" class="px-6 py-8 text-center">
                    <div class="flex flex-col items-center text-gray-500">
                        <i class="fas fa-microscope text-3xl mb-2"></i>
                        <p class="font-medium">No significant matches found</p>
                        <p class="text-sm mt-1">Try adjusting your search parameters</p>
                    </div>
                </td>
            </tr>
        `;
    }

     if (result.paths && result.paths.length > 0) {
    displayPathogenPaths(result.paths);
  } else {
    displayPathogenPaths([]);
  }


    setTimeout(() => {
        resultsSection.classList.remove('hidden');
        resultsSection.classList.add('animate-fade-in');
        resultsSection.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    }, 100);
}

// Helper function to extract sequence from GDN data
function extractSequenceForPosition(position, rawData) {
    if (!rawData) return "Sequence not available";
    
    try {

        const positionKey = position.replace(/:/g, '\\:'); // Escape special chars
        const regex = new RegExp(`${positionKey}\\|([A-Z]+)`);
        const match = rawData.match(regex);
        return match ? match[1] : "Sequence not found";
    } catch (e) {
        console.error("Error extracting sequence:", e);
        return "Error loading sequence";
    }
}

// score visualization 
function getScoreColor(score, type) {
    if (score >= 0.98) {
        return type === 'background' ? 'rgba(16, 185, 129, 0.7)' :
               type === 'border' ? 'rgba(16, 185, 129, 1)' : 'bg-green-500';
    }
    if (score >= 0.9) {
        return type === 'background' ? 'rgba(101, 163, 13, 0.7)' :
               type === 'border' ? 'rgba(101, 163, 13, 1)' : 'bg-lime-500';
    }
    if (score >= 0.8) {
        return type === 'background' ? 'rgba(234, 179, 8, 0.7)' :
               type === 'border' ? 'rgba(234, 179, 8, 1)' : 'bg-yellow-500';
    }
    if (score >= 0.6) {
        return type === 'background' ? 'rgba(249, 115, 22, 0.7)' :
               type === 'border' ? 'rgba(249, 115, 22, 1)' : 'bg-orange-500';
    }
    return type === 'background' ? 'rgba(239, 68, 68, 0.7)' :
           type === 'border' ? 'rgba(239, 68, 68, 1)' : 'bg-red-500';
}

function getScoreColorClass(score) {
    if (score >= 0.98) return 'bg-green-100 text-green-800';
    if (score >= 0.9) return 'bg-lime-100 text-lime-800';
    if (score >= 0.8) return 'bg-yellow-100 text-yellow-800';
    if (score >= 0.6) return 'bg-orange-100 text-orange-800';
    return 'bg-red-100 text-red-800';
}

function getScoreTextClass(score) {
    if (score >= 0.98) return 'text-green-700';
    if (score >= 0.9) return 'text-lime-700';
    if (score >= 0.8) return 'text-yellow-700';
    if (score >= 0.6) return 'text-orange-700';
    return 'text-red-700';
}

// Risk assessment logic
function determineRiskLevel(matches) {
    if (!matches.length) return 'None';
    if (matches.some(m => m.score >= 0.98)) return 'Critical';
    if (matches.some(m => m.score >= 0.9)) return 'High';
    if (matches.some(m => m.score >= 0.8)) return 'Moderate';
    if (matches.some(m => m.score >= 0.6)) return 'Low';
    return 'Minimal';
}

function getRiskLevelClass(matches) {
    const level = determineRiskLevel(matches);
    switch (level) {
        case 'Critical': return 'text-red-600';
        case 'High': return 'text-orange-500';
        case 'Moderate': return 'text-yellow-500';
        case 'Low': return 'text-lime-500';
        case 'Minimal': return 'text-green-400';
        default: return 'text-gray-500';
    }
}

// Global function to show match details
window.showMatchDetails = function(position, sequence) {
    const modal = document.createElement('div');
    modal.className = 'fixed inset-0 bg-black bg-opacity-75 flex items-center justify-center z-50 p-4';
    modal.innerHTML = `
        <div class="bg-white rounded-xl shadow-2xl max-w-4xl w-full max-h-[90vh] overflow-y-auto">
            <div class="sticky top-0 bg-white p-4 border-b flex justify-between items-center">
                <h3 class="text-xl font-bold text-gray-900">
                    <i class="fas fa-dna text-indigo-500 mr-2"></i>
                    Match Details: ${position}
                </h3>
                <button onclick="this.closest('div[style*=\"fixed\"]').remove()" 
                    class="text-gray-500 hover:text-gray-700">
                    <i class="fas fa-times text-xl"></i>
                </button>
            </div>
            
            <div class="p-6 space-y-6">
                <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
                    <div>
                        <h4 class="text-sm font-medium text-gray-500 mb-2">Sequence Information</h4>
                        <div class="bg-gray-50 p-4 rounded-lg">
                            <div class="font-mono text-sm overflow-x-auto">${sequence}</div>
                        </div>
                    </div>
                    
                    <div>
                        <h4 class="text-sm font-medium text-gray-500 mb-2">Match Characteristics</h4>
                        <div class="space-y-3">
                            <div class="flex justify-between">
                                <span class="text-gray-600">Position</span>
                                <span class="font-mono">${position}</span>
                            </div>
                            <div class="flex justify-between">
                                <span class="text-gray-600">Length</span>
                                <span>${sequence?.length || 'N/A'} bp</span>
                            </div>
                            <div class="flex justify-between">
                                <span class="text-gray-600">GC Content</span>
                                <span>${calculateGCContent(sequence)}%</span>
                            </div>
                        </div>
                    </div>
                </div>
                
                <div>
                    <h4 class="text-sm font-medium text-gray-500 mb-2">Sequence Alignment</h4>
                    <div class="bg-gray-900 text-green-400 p-4 rounded-lg font-mono text-xs overflow-x-auto">
                        <div class="mb-1">Reference: ${sequence}</div>
                        <div class="text-gray-400">${'|'.repeat(sequence?.length || 0)}</div>
                        <div>Query:     ${sequence}</div>
                    </div>
                </div>
            </div>
            
            <div class="sticky bottom-0 bg-gray-50 px-4 py-3 border-t flex justify-end space-x-3">
                <button onclick="this.closest('div[style*=\"fixed\"]').remove()" 
                    class="px-4 py-2 border border-gray-300 rounded-md text-gray-700 hover:bg-gray-100">
                    Close
                </button>
                <button onclick="exportMatch('${position}', '${sequence}')" 
                    class="px-4 py-2 bg-indigo-600 text-white rounded-md hover:bg-indigo-700">
                    <i class="fas fa-download mr-2"></i> Export
                </button>
            </div>
        </div>
    `;
    document.body.appendChild(modal);
};

// function to calculate gc function
function calculateGCContent(sequence) {
    if (!sequence) return 0;
    const gcCount = (sequence.match(/[GC]/gi) || []).length;
    return ((gcCount / sequence.length) * 100).toFixed(1);
}


  // Helper function for risk level colors
  function getRiskColorClass(riskLevel) {
    switch ((riskLevel || '').toLowerCase()) {
      case 'high': return 'text-red-600';
      case 'medium': return 'text-yellow-600';
      case 'low': return 'text-green-600';
      default: return 'text-gray-600';
    }
  }

  // Show error notification
  function showError(message) {
    const notification = document.createElement('div');
    notification.className = 'fixed top-4 right-4 bg-red-500 text-white px-4 py-2 rounded-lg shadow-lg flex items-center';
    notification.innerHTML = `
      <i class="fas fa-exclamation-circle mr-2"></i>
      <span>${message}</span>
      <button class="ml-4" onclick="this.parentElement.remove()">
        <i class="fas fa-times"></i>
      </button>
    `;
    document.body.appendChild(notification);
    
    // Auto-remove after 5 seconds
    setTimeout(() => notification.remove(), 5000);
  }
});