document.addEventListener('DOMContentLoaded', function() {
  const form = document.getElementById('similarityForm');
  const resultsSection = document.getElementById('resultsSection');
  const loadingSpinner = document.getElementById('loadingSpinner');
  const submitBtn = document.getElementById('submitBtn');
  const resetBtn = document.getElementById('resetBtn');
  
  // Form submission handler
  form.addEventListener('submit', async function(e) {
    e.preventDefault();
    
    try {
      // Show loading state
      submitBtn.disabled = true;
      submitBtn.innerHTML = '<i class="fas fa-spinner fa-spin mr-2"></i> Analyzing...';
      loadingSpinner.classList.remove('hidden');
      resultsSection.classList.add('hidden');
      
      // Get form data
      const formData = new FormData(e.target);
      const data = {
        dna1: formData.get('dna1'),
        dna2: formData.get('dna2'),
        kmerSize: parseInt(formData.get('kmerSize')),
        threshold: parseFloat(formData.get('threshold')),
        algorithm: formData.get('algorithm')
      };
      
      // Call backend API
      const response = await fetch("http://localhost:18080/calculate-similarity", {
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
      displayResults(result, data);

    } catch (error) {
      console.error("Analysis failed:", error);
      showError(error.message || "Analysis failed. Please try again.");
    } finally {
      // Reset loading state
      submitBtn.disabled = false;
      submitBtn.innerHTML = '<i class="fas fa-play mr-2"></i> Compare Sequences';
      loadingSpinner.classList.add('hidden');
    }
  });
  
  // Reset button handler
  resetBtn.addEventListener('click', function() {
    form.reset();
    resultsSection.classList.add('hidden');
  });
  
  // Display results function
  function displayResults(result, data) {
    const similarity = result.similarity;
    
    // Update main similarity score
    const similarityValue = document.getElementById('similarityValue');
    const similarityFill = document.getElementById('similarityFill');
    const similarityText = document.getElementById('similarityText');
    
    similarityValue.textContent = similarity.toFixed(3);
    similarityFill.style.width = `${similarity * 100}%`;
    
    // Update similarity text
    if (similarity >= 0.9) {
      similarityText.textContent = "Very High Similarity";
      similarityText.className = "mt-1 text-sm font-medium text-green-600";
    } else if (similarity >= 0.7) {
      similarityText.textContent = "High Similarity";
      similarityText.className = "mt-1 text-sm font-medium text-green-500";
    } else if (similarity >= 0.5) {
      similarityText.textContent = "Moderate Similarity";
      similarityText.className = "mt-1 text-sm font-medium text-yellow-600";
    } else {
      similarityText.textContent = "Low Similarity";
      similarityText.className = "mt-1 text-sm font-medium text-red-600";
    }
    
    // Update sequence stats
    document.getElementById('seqALength').textContent = data.dna1.length;
    document.getElementById('seqBLength').textContent = data.dna2.length;
    document.getElementById('seqAGC').textContent = calculateGCContent(data.dna1) + '%';
    document.getElementById('seqBGC').textContent = calculateGCContent(data.dna2) + '%';
    
    // Update algorithm
    document.getElementById('algorithmUsed').textContent = 
      data.algorithm === 'alignment' ? 'Sequence Alignment' : 'K-mer Similarity';
    
    // Update scores
    document.getElementById('kmerScore').textContent = result.kmerScore.toFixed(3);
    document.getElementById('kmerScoreBar').style.width = `${result.kmerScore * 100}%`;
    document.getElementById('alignmentScore').textContent = result.alignmentScore.toFixed(3);
    document.getElementById('alignmentScoreBar').style.width = `${result.alignmentScore * 100}%`;
    
    // Update match status
    const matchStatus = document.getElementById('matchStatus');
    if (result.match) {
      matchStatus.textContent = "SIGNIFICANT MATCH FOUND";
      matchStatus.className = "font-medium px-3 py-1 rounded-full bg-green-100 text-green-800 text-sm";
    } else {
      matchStatus.textContent = "NO SIGNIFICANT MATCH";
      matchStatus.className = "font-medium px-3 py-1 rounded-full bg-red-100 text-red-800 text-sm";
    }
    
    // Show results
    resultsSection.classList.remove('hidden');
    resultsSection.scrollIntoView({ behavior: 'smooth' });
  }
  
  // Helper function to calculate GC content
  function calculateGCContent(sequence) {
    if (!sequence) return 0;
    const gcCount = (sequence.match(/[GC]/gi) || []).length;
    return ((gcCount / sequence.length) * 100).toFixed(1);
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