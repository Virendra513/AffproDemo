<script>
  // Function to update the active class on sidebar items
  document.addEventListener("DOMContentLoaded", function() {
    // Get the active content ID from the hidden field
    const activeContentId = document.querySelector('input[name="activeContentId"]').value;

    // Set the corresponding sidebar item as active
    const sidebarItems = document.querySelectorAll('.unique-item');
    sidebarItems.forEach(item => {
      if (item.getAttribute('data-target') === activeContentId) {
        item.classList.add('active');
      } else {
        item.classList.remove('active');
      }
    });

    // Optionally, set the active content class as well
    const contentItems = document.querySelectorAll('.unique-content-item');
    contentItems.forEach(content => {
      if (content.id === 'content-4') {
        content.classList.add('active');
      } else {
        content.classList.remove('active');
      }
    });
  });
</script>
